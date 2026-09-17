function [Detail,Summary]=Run_BinnedTemplateNoiseAttribution(snrList,numSeeds,numLowRevs,numHighRevs,seedList)
%RUN_BINNEDTEMPLATENOISEATTRIBUTION Separate low/high/common noise effects.
% Uses the paper-baseline low-speed estimator: spatial-bin mean + SG smooth.
% The EO pair is fixed to isolate template calibration and parameter recovery.
if nargin<1||isempty(snrList),snrList=[10 5];end
if nargin<2||isempty(numSeeds),numSeeds=3;end
if nargin<3||isempty(numLowRevs),numLowRevs=20;end
if nargin<4||isempty(numHighRevs),numHighRevs=3;end
if nargin<5||isempty(seedList),seedList=1:numSeeds;else,numSeeds=numel(seedList);end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=numLowRevs;cfg0.NumRevs_high=numHighRevs;cfg0.route30ForwardModel="low_increment";cfg0.A_true=[.25 .15];fRot=cfg0.RPM_high/60;eoTrue=[10 26];cfg0.f_true=eoTrue*fRot;cfg0.phi_true=[pi/4 -pi/3];cfg0.route30DualFrequencyRangeHz=eoTrue*fRot;cfg0.dualSyncCandidateEOPairs=eoTrue;cfg0.route30GapHalfWidthMm=.20;cfg0.dualSyncGapCount=7;cfg0.dualSyncKeepPerGap=10;cfg0.dualSyncIteratedReplayCount=20;cfg0.dualSyncRefineCount=5;cfg0.structuredComponentAmplitudeFloorMm=.05;
g0=.50;delta=.10;truthT=repmat(eval_gap_template(lib,g0,lib.xGrid),1,numel(cfg0.alpha_k));truthCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg0.alpha_k));truthModel=build_path_template_model(lib,truthT,truthCal);
caseNames={'low_only','high_only','common'};nRows=numel(snrList)*numSeeds*numel(caseNames);Detail=repmat(empty_row(),nRows,1);ir=0;outDir=fullfile(root,'output','binned_template_noise_attribution');if ~exist(outDir,'dir'),mkdir(outDir);end;tag=sprintf('%dlevels_%dseeds_%dLrevs_%dHrevs',numel(snrList),numSeeds,numLowRevs,numHighRevs);
for isnr=1:numel(snrList)
    targetSnr=snrList(isnr);
    for icase=1:numel(caseNames)
        switch caseNames{icase}
            case 'low_only', lowSnr=targetSnr;highSnr=Inf;
            case 'high_only',lowSnr=Inf;highSnr=targetSnr;
            otherwise,lowSnr=targetSnr;highSnr=targetSnr;
        end
        for iseed=seedList
            % Pair the realizations across cases: low_only/common share the
            % same low-speed noise, while high_only/common share the same
            % high-speed noise. This makes source attribution interpretable.
            lowSeed=951000+1000*isnr+iseed;highSeed=952000+1000*isnr+iseed;
            cfg=cfg0;lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,lowSnr,lowSeed);high=simulate_highspeed_from_low_increment(truthModel,delta,cfg,highSnr,'snr_db',highSeed);highMap=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
            low=build_low_speed_templates_binned(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,'mean',struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',9,'minCoverage',.90));opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);[pc,cal]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;model=build_path_template_model(lib,low.templateBySensor,pc);state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfg);state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;fit=run_dual_sync_voltage_vp(highMap,model,cfg,state);
            ir=ir+1;Detail(ir)=pack_row(isnr,targetSnr,iseed,caseNames{icase},lowSnr,highSnr,low,truthT,lib.xGrid,pc,fit,g0,delta,cfg.A_true,cal);
        end
    end
    checkpoint=struct2table(Detail(1:ir));writetable(checkpoint,fullfile(outDir,['checkpoint_',tag,'.csv']));save(fullfile(outDir,['checkpoint_',tag,'.mat']),'Detail','ir','snrList','caseNames');
end
T=struct2table(Detail);Summary=make_summary(T,snrList,caseNames);writetable(T,fullfile(outDir,['detail_',tag,'.csv']));writetable(Summary,fullfile(outDir,['summary_',tag,'.csv']));save(fullfile(outDir,['attribution_',tag,'.mat']),'Detail','Summary','snrList','caseNames');plot_summary(Summary,snrList,caseNames,outDir);disp(Summary);
end

function r=empty_row()
r=struct('snr_index',NaN,'target_snr_db',NaN,'seed',NaN,'noise_case','','low_snr_db',NaN,'high_snr_db',NaN,'template_rmse_V',NaN,'weighted_derivative_rmse_V_per_mm',NaN,'derivative_correlation',NaN,'g0_error_mm',NaN,'delta_gap_error_mm',NaN,'A1_est_mm',NaN,'A2_est_mm',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,'Amax_error_mm',NaN,'success',false,'rmse_V',NaN,'calibration_rmse_V',NaN);
end
function r=pack_row(isnr,targetSnr,seed,caseName,lowSnr,highSnr,low,Ttrue,x,pc,fit,g0,delta,Atrue,cal)
d=median(diff(x));D=spatial_gradient(low.templateBySensor,d);Dt=spatial_gradient(Ttrue,d);de=D-Dt;w=Dt.^2;ae=fit.A_id(:);r=empty_row();r.snr_index=isnr;r.target_snr_db=targetSnr;r.seed=seed;r.noise_case=caseName;r.low_snr_db=lowSnr;r.high_snr_db=highSnr;r.template_rmse_V=sqrt(mean((low.templateBySensor-Ttrue).^2,'all'));r.weighted_derivative_rmse_V_per_mm=sqrt(sum(w.*de.^2,'all')/sum(w,'all'));r.derivative_correlation=sum(D.*Dt,'all')/sqrt(sum(D.^2,'all')*sum(Dt.^2,'all'));r.g0_error_mm=abs(pc.g0-g0);r.delta_gap_error_mm=abs((fit.g_used-pc.g0)-delta);r.A1_est_mm=ae(1);r.A2_est_mm=ae(2);r.A1_error_mm=abs(ae(1)-Atrue(1));r.A2_error_mm=abs(ae(2)-Atrue(2));r.Amax_error_mm=max(r.A1_error_mm,r.A2_error_mm);r.success=r.delta_gap_error_mm<=.01&&r.Amax_error_mm<=.01;r.rmse_V=fit.rmse;r.calibration_rmse_V=cal.rmse_V;
end
function S=make_summary(T,snrList,caseNames)
vars={'template_rmse_V','weighted_derivative_rmse_V_per_mm','derivative_correlation','g0_error_mm','delta_gap_error_mm','A1_error_mm','A2_error_mm','Amax_error_mm'};nRows=numel(snrList)*numel(caseNames);names=[{'snr_index','target_snr_db','noise_case','n'},reshape([strcat(vars,'_median');strcat(vars,'_iqr')],1,[]),{'success_rate'}];types=[{'double','double','string','double'},repmat({'double'},1,2*numel(vars)+1)];S=table('Size',[nRows,numel(names)],'VariableTypes',types,'VariableNames',names);row=0;
for i=1:numel(snrList),for c=1:numel(caseNames),row=row+1;idx=T.snr_index==i&strcmp(T.noise_case,caseNames{c});S.snr_index(row)=i;S.target_snr_db(row)=snrList(i);S.noise_case(row)=string(caseNames{c});S.n(row)=nnz(idx);for iv=1:numel(vars),v=vars{iv};q=T.(v)(idx);S.([v,'_median'])(row)=median(q);S.([v,'_iqr'])(row)=iqr(q);end;S.success_rate(row)=mean(T.success(idx));end,end
end
function plot_summary(S,snrList,caseNames,outDir)
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 13]);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');cc=lines(numel(caseNames));labels=arrayfun(@(x)sprintf('%g',x),snrList,'UniformOutput',false);panels={{'weighted_derivative_rmse_V_per_mm_median','weighted_derivative_rmse_V_per_mm_iqr',1,'Weighted derivative RMSE (V/mm)'},{'delta_gap_error_mm_median','delta_gap_error_mm_iqr',1000,'Clearance-change error (\mum)'},{'Amax_error_mm_median','Amax_error_mm_iqr',1000,'Maximum amplitude error (\mum)'},{'success_rate','',100,'Success rate (%)'}};
for ip=1:4,nexttile;hold on;box on;for c=1:numel(caseNames),idx=strcmp(S.noise_case,caseNames{c});[~,o]=sort(S.snr_index(idx));q=find(idx);q=q(o);y=panels{ip}{3}*S.(panels{ip}{1})(q);if isempty(panels{ip}{2}),plot(1:numel(snrList),y,'-o','Color',cc(c,:),'LineWidth',1.1,'MarkerSize',4,'DisplayName',strrep(caseNames{c},'_',' '));else,errorbar(1:numel(snrList),y,.5*panels{ip}{3}*S.(panels{ip}{2})(q),'-o','Color',cc(c,:),'LineWidth',1.1,'MarkerSize',4,'DisplayName',strrep(caseNames{c},'_',' '));end;end;set(gca,'XTick',1:numel(snrList),'XTickLabel',labels);xlabel('Target SNR (dB)');ylabel(panels{ip}{4});if ip==1,legend('Location','best');end;end
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');set(findall(fig,'-property','FontSize'),'FontSize',8);set(findall(fig,'Type','axes'),'TickDir','in');exportgraphics(fig,fullfile(outDir,'noise_attribution.png'),'Resolution',240);try,exportgraphics(fig,fullfile(outDir,'noise_attribution.pdf'),'ContentType','vector');catch,end
end
function D=spatial_gradient(T,d)
if isvector(T),D=gradient(T,d);else,[~,D]=gradient(T,d,d);end
end
