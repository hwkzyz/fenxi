function [Detail,Summary]=Run_BinnedTemplateCommonSNRSweep(snrList,numSeeds,numLowRevs,numHighRevs,methodList)
%RUN_BINNEDTEMPLATECOMMONSNRSWEEP Common low/high SNR sweep.
% Low- and high-speed records use the same SNR value but independent noise.
% EO is fixed here to isolate template/parameter estimation; full EO search
% remains a separate validation step.
if nargin<1||isempty(snrList),snrList=[Inf 30 20 15 10 5];end
if nargin<2||isempty(numSeeds),numSeeds=3;end
if nargin<3||isempty(numLowRevs),numLowRevs=20;end
if nargin<4||isempty(numHighRevs),numHighRevs=8;end
if nargin<5||isempty(methodList),methodList={'mean','median','huber'};end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=numLowRevs;cfg0.NumRevs_high=numHighRevs;cfg0.route30ForwardModel="low_increment";cfg0.A_true=[.25 .15];fRot=cfg0.RPM_high/60;eoTrue=[10 26];cfg0.f_true=eoTrue*fRot;cfg0.phi_true=[pi/4 -pi/3];cfg0.route30DualFrequencyRangeHz=eoTrue*fRot;cfg0.dualSyncCandidateEOPairs=eoTrue;cfg0.route30GapHalfWidthMm=.20;cfg0.dualSyncGapCount=7;cfg0.dualSyncKeepPerGap=10;cfg0.dualSyncIteratedReplayCount=20;cfg0.dualSyncRefineCount=5;cfg0.structuredComponentAmplitudeFloorMm=.05;
g0=.50;delta=.10;truthT=repmat(eval_gap_template(lib,g0,lib.xGrid),1,numel(cfg0.alpha_k));truthCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg0.alpha_k));truthModel=build_path_template_model(lib,truthT,truthCal);methods=cellstr(methodList);
Detail=repmat(empty_row(),numel(snrList)*numSeeds*numel(methods),1);ir=0;outDir=fullfile(root,'output','binned_template_common_snr_sweep');if ~exist(outDir,'dir'),mkdir(outDir);end;tag=sprintf('%dlevels_%dseeds_%dLrevs_%dHrevs_%dmethods',numel(snrList),numSeeds,numLowRevs,numHighRevs,numel(methods));
for isnr=1:numel(snrList)
    snrDb=snrList(isnr);
    for iseed=1:numSeeds
        cfg=cfg0;lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,snrDb,941000+100*isnr+iseed);high=simulate_highspeed_from_low_increment(truthModel,delta,cfg,snrDb,'snr_db',942000+100*isnr+iseed);highMap=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
        for im=1:numel(methods)
            low=build_low_speed_templates_binned(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,methods{im},struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',9,'minCoverage',.90));opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,cal]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;model=build_path_template_model(lib,low.templateBySensor,pc);state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;fit=run_dual_sync_voltage_vp(highMap,model,cfg,state);ir=ir+1;Detail(ir)=pack_row(isnr,snrDb,iseed,methods{im},low,truthT,lib.xGrid,pc,fit,g0,delta,cfg.A_true,cal);
        end
        checkpoint=struct2table(Detail(1:ir));writetable(checkpoint,fullfile(outDir,['checkpoint_',tag,'.csv']));save(fullfile(outDir,['checkpoint_',tag,'.mat']),'Detail','ir','snrList','methods');
    end
end
T=struct2table(Detail);Summary=make_summary(T,snrList,methods);writetable(T,fullfile(outDir,['detail_',tag,'.csv']));writetable(Summary,fullfile(outDir,['summary_',tag,'.csv']));save(fullfile(outDir,['sweep_',tag,'.mat']),'Detail','Summary','snrList','methods');plot_summary(Summary,snrList,methods,outDir);disp(Summary);
end

function r=empty_row()
r=struct('snr_index',NaN,'snr_db',NaN,'seed',NaN,'method','','template_rmse_V',NaN,'weighted_derivative_rmse_V_per_mm',NaN,'derivative_correlation',NaN,'g0_error_mm',NaN,'delta_gap_error_mm',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,'Amax_error_mm',NaN,'success',false,'rmse_V',NaN,'calibration_rmse_V',NaN);
end
function r=pack_row(isnr,snrDb,seed,name,low,Ttrue,x,pc,fit,g0,delta,Atrue,cal)
d=median(diff(x));D=spatial_gradient(low.templateBySensor,d);Dt=spatial_gradient(Ttrue,d);de=D-Dt;w=Dt.^2;ae=fit.A_id(:);r=empty_row();r.snr_index=isnr;r.snr_db=snrDb;r.seed=seed;r.method=name;r.template_rmse_V=sqrt(mean((low.templateBySensor-Ttrue).^2,'all'));r.weighted_derivative_rmse_V_per_mm=sqrt(sum(w.*de.^2,'all')/sum(w,'all'));r.derivative_correlation=sum(D.*Dt,'all')/sqrt(sum(D.^2,'all')*sum(Dt.^2,'all'));r.g0_error_mm=abs(pc.g0-g0);r.delta_gap_error_mm=abs((fit.g_used-pc.g0)-delta);r.A1_error_mm=abs(ae(1)-Atrue(1));r.A2_error_mm=abs(ae(2)-Atrue(2));r.Amax_error_mm=max(r.A1_error_mm,r.A2_error_mm);r.success=r.delta_gap_error_mm<=.01&&r.Amax_error_mm<=.01;r.rmse_V=fit.rmse;r.calibration_rmse_V=cal.rmse_V;
end
function S=make_summary(T,snrList,methods)
vars={'template_rmse_V','weighted_derivative_rmse_V_per_mm','derivative_correlation','g0_error_mm','delta_gap_error_mm','A1_error_mm','A2_error_mm','Amax_error_mm'};
nRows=numel(snrList)*numel(methods);S=table('Size',[nRows,4+2*numel(vars)+1], ...
    'VariableTypes',[{'double','double','string','double'},repmat({'double'},1,2*numel(vars)+1)], ...
    'VariableNames',[{'snr_index','snr_db','method','n'},reshape([strcat(vars,'_median');strcat(vars,'_iqr')],1,[]),{'success_rate'}]);row=0;
for i=1:numel(snrList),for m=1:numel(methods),row=row+1;idx=T.snr_index==i&strcmp(T.method,methods{m});S.snr_index(row)=i;S.snr_db(row)=snrList(i);S.method(row)=string(methods{m});S.n(row)=nnz(idx);for iv=1:numel(vars),v=vars{iv};q=T.(v)(idx);S.([v,'_median'])(row)=median(q);S.([v,'_iqr'])(row)=iqr(q);end;S.success_rate(row)=mean(T.success(idx));end,end
end
function plot_summary(S,snrList,methods,outDir)
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 13]);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');cc=lines(numel(methods));labels=arrayfun(@snr_label,snrList,'UniformOutput',false);panels={{'template_rmse_V_median','template_rmse_V_iqr',1000,'Template RMSE (mV)'},{'weighted_derivative_rmse_V_per_mm_median','weighted_derivative_rmse_V_per_mm_iqr',1,'Weighted derivative RMSE (V/mm)'},{'delta_gap_error_mm_median','delta_gap_error_mm_iqr',1000,'Clearance-change error (\mum)'},{'Amax_error_mm_median','Amax_error_mm_iqr',1000,'Maximum amplitude error (\mum)'}};
for ip=1:4,nexttile;hold on;box on;for m=1:numel(methods),idx=strcmp(S.method,methods{m});[~,o]=sort(S.snr_index(idx));q=find(idx);q=q(o);errorbar(1:numel(snrList),panels{ip}{3}*S.(panels{ip}{1})(q),.5*panels{ip}{3}*S.(panels{ip}{2})(q),'-o','Color',cc(m,:),'LineWidth',1.1,'MarkerSize',4,'DisplayName',methods{m});end;set(gca,'XTick',1:numel(snrList),'XTickLabel',labels);xlabel('Common low/high SNR (dB)');ylabel(panels{ip}{4});if ip==1,legend('Location','best');end;end
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');set(findall(fig,'-property','FontSize'),'FontSize',8);set(findall(fig,'Type','axes'),'TickDir','in');exportgraphics(fig,fullfile(outDir,'common_snr_sweep.png'),'Resolution',240);try,exportgraphics(fig,fullfile(outDir,'common_snr_sweep.pdf'),'ContentType','vector');end
end
function s=snr_label(x)
if isinf(x),s='Noiseless';else,s=sprintf('%g',x);end
end
function D=spatial_gradient(T,d)
if isvector(T),D=gradient(T,d);else,[~,D]=gradient(T,d,d);end
end
