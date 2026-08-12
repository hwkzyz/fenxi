function [Candidate,NoiseSource,Spectrum]=Run_EOCombinationOrderDiagnostic()
%RUN_EOCOMBINATIONORDERDIAGNOSTIC Diagnose 16/20 order competition.
% Tests whether nonlinear combination orders exist without noise and whether
% low/high-speed noise makes their candidate basins overtake EO=(10,26).
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;cfg0.route30ForwardModel="low_increment";cfg0.A_true=[.25 .15];fRot=cfg0.RPM_high/60;eoTrue=[10 26];cfg0.f_true=eoTrue*fRot;cfg0.phi_true=[pi/4 -pi/3];cfg0.route30DualFrequencyRangeHz=eoTrue*fRot;cfg0.route30GapHalfWidthMm=.20;cfg0.dualSyncGapCount=11;cfg0.dualSyncKeepPerGap=80;cfg0.dualSyncIteratedReplayCount=200;cfg0.dualSyncRefineCount=20;cfg0.structuredComponentAmplitudeFloorMm=.05;
g0=.50;delta=.10;truthT=repmat(eval_gap_template(lib,g0,lib.xGrid),1,numel(cfg0.alpha_k));truthCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg0.alpha_k));truthModel=build_path_template_model(lib,truthT,truthCal);
snrList=[Inf 25 15 5];pairs=[10 26;16 26;16 20;20 22];rows=repmat(empty_candidate_row(),numel(snrList)*size(pairs,1),1);ir=0;noiselessMap=[];
for is=1:numel(snrList)
    snrDb=snrList(is);[model,state,highMap,low]=make_case(snrDb,snrDb,941001,942001,cfg0,lib,truthModel,g0,delta);
    if isinf(snrDb),noiselessMap=highMap;end
    for ip=1:size(pairs,1)
        cfg=cfg0;cfg.dualSyncCandidateEOPairs=pairs(ip,:);fit=run_dual_sync_voltage_vp(highMap,model,cfg,state);ir=ir+1;rows(ir)=pack_candidate_row(snrDb,pairs(ip,:),fit,low,truthT,lib.xGrid,g0,delta,cfg.A_true);
    end
end
Candidate=struct2table(rows);Candidate.rmse_excess_V=zeros(height(Candidate),1);for is=1:numel(snrList),idx=find(Candidate.snr_index==is);Candidate.rmse_excess_V(idx)=Candidate.rmse_V(idx)-min(Candidate.rmse_V(idx));end

% Separate low-speed and high-speed 5 dB noise using paired realizations.
caseNames={'low_only','high_only','common'};src=repmat(empty_source_row(),numel(caseNames),1);
for ic=1:numel(caseNames)
    switch caseNames{ic},case 'low_only',ls=5;hs=Inf;case 'high_only',ls=Inf;hs=5;otherwise,ls=5;hs=5;end
    [model,state,highMap]=make_case(ls,hs,941001,942001,cfg0,lib,truthModel,g0,delta);cfg=cfg0;cfg.dualSyncCandidateEOPairs=[];fit=run_dual_sync_voltage_vp(highMap,model,cfg,state);src(ic)=pack_source_row(caseNames{ic},fit,state,g0,delta,cfg.A_true,eoTrue);
end
NoiseSource=struct2table(src);

% Noiseless displacement-proxy spectrum. Division by local sensitivity
% removes most spatial amplitude modulation but retains nonlinear terms.
sid=noiselessMap.S_v(:);x=noiselessMap.x_v(:);theta=noiselessMap.theta_v(:);V=noiselessMap.V_a(:);V0=eval_gap_template(truthModel,g0+delta,x,sid);Fx=eval_gap_derivative(truthModel,g0+delta,x,sid);thr=prctile(abs(Fx(isfinite(Fx))),40);good=isfinite(V)&isfinite(V0)&isfinite(Fx)&abs(Fx)>=thr;uProxy=-(V(good)-V0(good))./Fx(good);th=theta(good);maxEO=40;X=zeros(numel(th),2*maxEO);for k=1:maxEO,X(:,2*k-1)=sin(k*th);X(:,2*k)=cos(k*th);end;coef=X\uProxy;amp=zeros(maxEO,1);for k=1:maxEO,amp(k)=hypot(coef(2*k-1),coef(2*k));end;Spectrum=table((1:maxEO).',amp,'VariableNames',{'EO','proxy_amplitude_mm'});

outDir=fullfile(root,'output','eo_combination_order_diagnostic');if ~exist(outDir,'dir'),mkdir(outDir);end;writetable(Candidate,fullfile(outDir,'candidate_pair_comparison.csv'));writetable(NoiseSource,fullfile(outDir,'noise_source_full_eo.csv'));writetable(Spectrum,fullfile(outDir,'noiseless_proxy_order_spectrum.csv'));save(fullfile(outDir,'diagnostic.mat'),'Candidate','NoiseSource','Spectrum','pairs','snrList');plot_results(Candidate,Spectrum,snrList,pairs,outDir);disp(Candidate);disp(NoiseSource);disp(Spectrum(ismember(Spectrum.EO,[10 16 20 26]),:));
end

function [model,state,highMap,low]=make_case(lowSnr,highSnr,lowSeed,highSeed,cfg,lib,truthModel,g0,delta)
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,lowSnr,lowSeed);high=simulate_highspeed_from_low_increment(truthModel,delta,cfg,highSnr,'snr_db',highSeed);highMap=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);low=build_low_speed_templates_binned(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,'mean',struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',9,'minCoverage',.90));opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,~]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;model=build_path_template_model(lib,low.templateBySensor,pc);state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;
end
function r=empty_candidate_row()
r=struct('snr_index',NaN,'snr_db',NaN,'eo1',NaN,'eo2',NaN,'rmse_V',NaN,'rmse_excess_V',NaN,'delta_gap_est_mm',NaN,'delta_gap_error_mm',NaN,'A1_est_mm',NaN,'A2_est_mm',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,'template_derivative_correlation',NaN);
end
function r=pack_candidate_row(snrDb,pair,fit,low,Ttrue,x,g0,delta,Atrue)
d=median(diff(x));D=spatial_gradient(low.templateBySensor,d);Dt=spatial_gradient(Ttrue,d);r=empty_candidate_row();if isinf(snrDb),r.snr_index=1;elseif snrDb==25,r.snr_index=2;elseif snrDb==15,r.snr_index=3;else,r.snr_index=4;end;r.snr_db=snrDb;r.eo1=pair(1);r.eo2=pair(2);r.rmse_V=fit.rmse;r.delta_gap_est_mm=fit.g_used-(g0+(fit.staticState.gHat-g0));r.delta_gap_error_mm=abs((fit.g_used-fit.staticState.gHat)-delta);r.A1_est_mm=fit.A_id(1);r.A2_est_mm=fit.A_id(2);r.A1_error_mm=abs(fit.A_id(1)-Atrue(1));r.A2_error_mm=abs(fit.A_id(2)-Atrue(2));r.template_derivative_correlation=sum(D.*Dt,'all')/sqrt(sum(D.^2,'all')*sum(Dt.^2,'all'));
end
function r=empty_source_row()
r=struct('noise_case','','eo1_est',NaN,'eo2_est',NaN,'eo_error',NaN,'A1_est_mm',NaN,'A2_est_mm',NaN,'delta_gap_est_mm',NaN,'delta_gap_error_mm',NaN,'rmse_V',NaN,'identification_confident',false,'noise_normalized_margin',NaN);
end
function r=pack_source_row(name,fit,state,g0,delta,Atrue,eoTrue)
r=empty_source_row();r.noise_case=name;r.eo1_est=fit.eo_id(1);r.eo2_est=fit.eo_id(2);r.eo_error=max(abs(fit.eo_id(:)-eoTrue(:)));r.A1_est_mm=fit.A_id(1);r.A2_est_mm=fit.A_id(2);r.delta_gap_est_mm=fit.g_used-state.gHat;r.delta_gap_error_mm=abs(r.delta_gap_est_mm-delta);r.rmse_V=fit.rmse;r.identification_confident=fit.identification_confident;r.noise_normalized_margin=fit.noise_normalized_margin; %#ok<NASGU>
end
function plot_results(C,S,snrList,pairs,outDir)
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 8]);tiledlayout(1,2,'TileSpacing','compact','Padding','compact');nexttile;hold on;box on;cc=lines(size(pairs,1));for ip=1:size(pairs,1),idx=C.eo1==pairs(ip,1)&C.eo2==pairs(ip,2);plot(1:numel(snrList),1000*C.rmse_excess_V(idx),'-o','LineWidth',1.1,'MarkerSize',4,'Color',cc(ip,:),'DisplayName',sprintf('(%d,%d)',pairs(ip,1),pairs(ip,2)));end;set(gca,'XTick',1:numel(snrList),'XTickLabel',{'Noiseless','25','15','5'});xlabel('Common low/high SNR (dB)');ylabel('RMSE excess above best (mV)');legend('Location','best');nexttile;bar(S.EO,1000*S.proxy_amplitude_mm,'FaceColor',[.2 .45 .75]);hold on;for e=[10 16 20 26],xline(e,'--','LineWidth',.8);end;xlim([0 40]);xlabel('Engine order');ylabel('Noiseless proxy amplitude (\mum)');box on;set(findall(fig,'-property','FontName'),'FontName','Times New Roman');set(findall(fig,'-property','FontSize'),'FontSize',8);set(findall(fig,'Type','axes'),'TickDir','in');exportgraphics(fig,fullfile(outDir,'eo_combination_diagnostic.png'),'Resolution',240);try,exportgraphics(fig,fullfile(outDir,'eo_combination_diagnostic.pdf'),'ContentType','vector');catch,end
end
function D=spatial_gradient(T,d)
if isvector(T),D=gradient(T,d);else,[~,D]=gradient(T,d,d);end
end
