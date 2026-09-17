function Results = Run_LowSpeedTemplateMethodImpact(lowSnrDb,highSnrDb,numLowRevs)
%RUN_LOWSPEEDTEMPLATEMETHODIMPACT Quantify template-processing impact.
% Noise-only experiment: the same noisy low-speed record and the same noisy
% high-speed record are inverted with several low-speed template estimators.
if nargin<1||isempty(lowSnrDb),lowSnrDb=15;end
if nargin<2||isempty(highSnrDb),highSnrDb=15;end
if nargin<3||isempty(numLowRevs),numLowRevs=20;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
demo=make_inv_log_2_demo_case(ctx,.2);cfg=demo.cfgCase;cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=numLowRevs;cfg.NumRevs_high=8;
cfg.route30ForwardModel="low_increment";cfg.A_true=[.25 .15];fRot=cfg.RPM_high/60;eoTrue=[10 26];cfg.f_true=eoTrue*fRot;cfg.phi_true=[pi/4 -pi/3];cfg.route30DualFrequencyRangeHz=eoTrue*fRot;cfg.dualSyncCandidateEOPairs=eoTrue;cfg.dualSyncGapCount=5;cfg.route30GapHalfWidthMm=.20;cfg.dualSyncKeepPerGap=5;cfg.dualSyncIteratedReplayCount=5;cfg.dualSyncRefineCount=2;cfg.structuredComponentAmplitudeFloorMm=.05;
g0=.50;delta=.10;
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,lowSnrDb,941001);
truthT=eval_gap_template(lib,g0,lib.xGrid);truthT=repmat(truthT,1,numel(cfg.alpha_k));trueCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));truthModel=build_path_template_model(lib,truthT,trueCal);
high=simulate_highspeed_from_low_increment(truthModel,delta,cfg,highSnrDb,'snr_db',942001);highMap=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
names={'interp_mean_sg','interp_median_sg','bin_mean_sg','bin_median_sg','bin_huber_sg'};rows=repmat(empty_row(),numel(names),1);templateStack=nan(numel(lib.xGrid),numel(names));derivativeStack=templateStack;
for im=1:numel(names)
    if startsWith(names{im},'bin_')
        statistic=extractBetween(names{im},'bin_','_sg');low=build_low_speed_templates_binned(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,char(statistic),struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',9,'minCoverage',.90));
    else
        method=strrep(names{im},'interp_','');low=aggregate_low_speed_template_methods(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,method,struct('smoothSpan',9));
    end
    opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,cal]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;model=build_path_template_model(lib,low.templateBySensor,pc);
    state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;fit=run_dual_sync_voltage_vp(highMap,model,cfg,state);diag=oracle_diagnostic(highMap,model,truthModel,pc,cfg,delta);
    d=median(diff(lib.xGrid));rows(im)=pack_row(names{im},low,truthT,lib.xGrid,pc,fit,g0,g0+delta,cfg.A_true,eoTrue,cal,diag);
    templateStack(:,im)=low.templateBySensor(:,1);derivativeStack(:,im)=spatial_gradient(templateStack(:,im),d);
end
Results=struct2table(rows);outDir=fullfile(root,'output','low_speed_template_method_impact');if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Results,fullfile(outDir,sprintf('method_impact_low_%gdb_high_%gdb_%drevs.csv',lowSnrDb,highSnrDb,numLowRevs)));
save(fullfile(outDir,sprintf('method_impact_low_%gdb_high_%gdb_%drevs.mat',lowSnrDb,highSnrDb,numLowRevs)),'Results','templateStack','derivativeStack','truthT','lib','names');
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 10]);tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;bar(categorical(Results.method),1000*Results.template_rmse_V);ylabel('Template RMSE (mV)');title('Template voltage error');box on;
nexttile;bar(categorical(Results.method),Results.derivative_rmse_V_per_mm);ylabel('Derivative RMSE (V/mm)');title('Template derivative error');box on;
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');set(findall(fig,'-property','FontSize'),'FontSize',8);set(findall(fig,'Type','axes'),'TickDir','in');exportgraphics(fig,fullfile(outDir,'method_impact_summary.png'),'Resolution',240);try,exportgraphics(fig,fullfile(outDir,'method_impact_summary.pdf'),'ContentType','vector');end
fig2=figure('Color','w','Units','centimeters','Position',[2 2 17 12]);tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;hold on;box on;cc=lines(numel(names));for im=1:numel(names),plot(lib.xGrid,templateStack(:,im),'Color',cc(im,:),'LineWidth',1.0,'DisplayName',names{im});end;plot(lib.xGrid,truthT(:,1),'k--','LineWidth',1.2,'DisplayName','truth');xlabel('Spatial coordinate x (mm)');ylabel('Template voltage (V)');title('Low-speed template comparison');legend('Location','best');
nexttile;hold on;box on;for im=1:numel(names),plot(lib.xGrid,derivativeStack(:,im),'Color',cc(im,:),'LineWidth',1.0,'DisplayName',names{im});end;plot(lib.xGrid,spatial_gradient(truthT(:,1),median(diff(lib.xGrid))),'k--','LineWidth',1.2,'DisplayName','truth');xlabel('Spatial coordinate x (mm)');ylabel('dT/dx (V/mm)');title('Template derivative comparison');legend('Location','best');
set(findall(fig2,'-property','FontName'),'FontName','Times New Roman');set(findall(fig2,'-property','FontSize'),'FontSize',8);set(findall(fig2,'Type','axes'),'TickDir','in');exportgraphics(fig2,fullfile(outDir,'method_impact_waveforms.png'),'Resolution',240);try,exportgraphics(fig2,fullfile(outDir,'method_impact_waveforms.pdf'),'ContentType','vector');end
disp(Results);
end

function r=empty_row()
r=struct('method','','template_rmse_V',NaN,'derivative_rmse_V_per_mm',NaN,'weighted_derivative_rmse_V_per_mm',NaN,'derivative_correlation',NaN,'derivative_gain',NaN,'total_sensitivity_rmse_V_per_mm',NaN,'g0_error_mm',NaN,'gap_error_mm',NaN,'delta_gap_est_mm',NaN,'delta_gap_error_mm',NaN,'eo_est_1',NaN,'eo_est_2',NaN,'eo_error',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,'A1_est_mm',NaN,'A2_est_mm',NaN,'oracle_A1_est_mm',NaN,'oracle_A2_est_mm',NaN,'oracle_delta_gap_est_mm',NaN,'mse_true_parameters_V2',NaN,'mse_normal_opt_V2',NaN,'mse_oracle_opt_V2',NaN,'rmse_V',NaN,'calibration_rmse_V',NaN,'success',false);
end
function r=pack_row(name,low,Ttrue,x,pc,fit,g0,gHigh,Atrue,eotrue,cal,diag)
d=median(diff(x));D=spatial_gradient(low.templateBySensor,d);Dt=spatial_gradient(Ttrue,d);de=D-Dt;w=Dt.^2;e=sqrt(mean((low.templateBySensor-Ttrue).^2,'all'));ed=sqrt(mean(de.^2,'all'));ew=sqrt(sum(w.*de.^2,'all')/sum(w,'all'));rho=sum(D.*Dt,'all')/sqrt(sum(D.^2,'all')*sum(Dt.^2,'all'));gain=sum(D.*Dt,'all')/sum(Dt.^2,'all');eo=fit.eo_id(:);ae=fit.A_id(:);r=empty_row();r.method=name;r.template_rmse_V=e;r.derivative_rmse_V_per_mm=ed;r.weighted_derivative_rmse_V_per_mm=ew;r.derivative_correlation=rho;r.derivative_gain=gain;r.total_sensitivity_rmse_V_per_mm=diag.total_sensitivity_rmse;r.g0_error_mm=abs(pc.g0-g0);r.gap_error_mm=abs(fit.g_used-gHigh);r.delta_gap_est_mm=fit.g_used-pc.g0;r.delta_gap_error_mm=abs(r.delta_gap_est_mm-(gHigh-g0));r.eo_est_1=eo(1);r.eo_est_2=eo(2);r.eo_error=max(abs(eo-eotrue(:)));r.A1_error_mm=abs(ae(1)-Atrue(1));r.A2_error_mm=abs(ae(2)-Atrue(2));r.A1_est_mm=ae(1);r.A2_est_mm=ae(2);r.oracle_A1_est_mm=diag.oracle_A(1);r.oracle_A2_est_mm=diag.oracle_A(2);r.oracle_delta_gap_est_mm=diag.oracle_delta;r.mse_true_parameters_V2=diag.mse_true;r.mse_normal_opt_V2=fit.rmse^2;r.mse_oracle_opt_V2=diag.mse_oracle;r.rmse_V=fit.rmse;r.calibration_rmse_V=cal.rmse_V;r.success=r.eo_error==0&&r.delta_gap_error_mm<=.05&&max([r.A1_error_mm,r.A2_error_mm])<=.05;
end

function D=spatial_gradient(T,d)
if isvector(T)
    D=gradient(T,d);
else
    [~,D]=gradient(T,d,d);
end
end

function D=oracle_diagnostic(map,model,truthModel,pc,cfg,delta)
theta=map.theta_v(:);x=map.x_v(:);V=map.V_a(:);sid=map.S_v(:);eo=[10 26];A=cfg.A_true(:).';phi=cfg.phi_true(:).';coef=[A(1)*cos(phi(1)),A(1)*sin(phi(1)),A(2)*cos(phi(2)),A(2)*sin(phi(2))];z0=[pc.g0+delta,0,coef];
predTrue=model_voltage(z0,eo,theta,x,sid,model);rTrue=safe_residual(V,predTrue);D.mse_true=mean(rTrue.^2);
lb=z0+[-.05,-.15,-.12,-.12,-.12,-.12];ub=z0+[.05,.15,.12,.12,.12,.12];[z,~]=solve_lsq_bounded(@(q)safe_residual(V,model_voltage(q,eo,theta,x,sid,model)),z0,lb,ub,100,1e-12,1e-12);pred=model_voltage(z,eo,theta,x,sid,model);r=safe_residual(V,pred);D.mse_oracle=mean(r.^2);D.oracle_A=[hypot(z(3),z(4)),hypot(z(5),z(6))];D.oracle_delta=z(1)-pc.g0;
xg=model.xGrid(:);h=median(diff(xg));dest=gradient(eval_gap_template(model,pc.g0+delta,xg,1),h);dtrue=gradient(eval_gap_template(truthModel,.5+delta,xg,1),h);D.total_sensitivity_rmse=sqrt(mean((dest-dtrue).^2));
end
function y=model_voltage(z,eo,theta,x,sid,model)
u=z(3)*sin(eo(1)*theta)+z(4)*cos(eo(1)*theta)+z(5)*sin(eo(2)*theta)+z(6)*cos(eo(2)*theta);y=eval_gap_template(model,z(1),x-z(2)-u,sid);
end
function r=safe_residual(V,pred)
r=V-pred;bad=~isfinite(r);r(bad)=10*max(std(V),1e-3);
end
