function [Detail,Summary,ReferenceNoise]=Run_PaperBaselineFixedNoiseEOStrategy(numSeeds,snrList)
%RUN_PAPERBASELINEFIXEDNOISEEOSTRATEGY Unified paper-baseline pilot.
% Low and high records use the same fixed absolute voltage-noise standard
% deviation, but independent random samples.  A noisy low-speed record is
% calibrated once and reused by strategies A and C on exactly the same
% high-speed record.  The high-speed truth is generated from an independent
% noise-free low-speed truth template, not from the estimated noisy template.
if nargin<1||isempty(numSeeds),numSeeds=3;end
if nargin<2||isempty(snrList),snrList=[25 15 5];end

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','paper_baseline_fixed_noise_eo_strategy');
if ~exist(outDir,'dir'),mkdir(outDir);end

ctx=load_inv_log_2_project_context(mainDir);
tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);
cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="low_increment";
cfg0.route30LowTemplateMethod="adaptive_sg";
cfg0.route30LowTemplateOptions=struct('gridSpacingMm',.02,...
    'minBinCount',5,'minCoverage',.90,'numFolds',5,...
    'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg0.route30FrequencyStructureMode="dual_sync_sync";
cfg0.route30StructuredFastBudgets=false;
cfg0.A_true=[.25 .15];cfg0.phi_true=[pi/4 -pi/3];
fRot=cfg0.RPM_high/60;eoTrue=[10 26];cfg0.f_true=eoTrue*fRot;
cfg0.route30DualFrequencyRangeHz=[min(cfg0.f_true) max(cfg0.f_true)];
cfg0.route30GapHalfWidthMm=.20;cfg0.dualSyncGapCount=11;
cfg0.dualSyncKeepPerGap=80;cfg0.dualSyncIteratedReplayCount=2100;
cfg0.dualSyncRefineCount=50;cfg0.dualSyncTopPairCount=40;
cfg0.dualSyncStartsPerTopPair=2;cfg0.structuredComponentAmplitudeFloorMm=.05;
cfg0.returnDualSyncDiagnostics=true;

gLow=.50;deltaTrue=.10;gHigh=gLow+deltaTrue;
truthT=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg0.alpha_k));
truthCal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,...
    'sensorIds',1:numel(cfg0.alpha_k));
truthModel=build_path_template_model(lib,truthT,truthCal);
sigmaRef=reference_noise_std(truthModel,deltaTrue,cfg0,lib.domain,snrList);
ReferenceNoise=table(snrList(:),sigmaRef(:),...
    'VariableNames',{'reference_snr_db','fixed_noise_std_V'});
writetable(ReferenceNoise,fullfile(outDir,'reference_noise.csv'));

strategyKeys=["global_budget","coarse_top_m_multistart"];
strategyLabels=["A: global-50","C: Top40 x 2"];
nMax=numSeeds*numel(snrList)*numel(strategyKeys);
rows=repmat(empty_row(),nMax,1);ir=0;
checkpoint=fullfile(outDir,'detail_checkpoint.csv');
if exist(checkpoint,'file')
    old=readtable(checkpoint,'TextType','string');
    keep=ismember(old.seed,1:numSeeds)&ismember(old.reference_snr_db,snrList)&...
        ismember(old.strategy_key,strategyKeys);
    oldRows=table2struct(old(keep,:));nOld=min(numel(oldRows),nMax);
    rows(1:nOld)=oldRows(1:nOld);ir=nOld;
    fprintf('Resuming from %d completed runs.\n',ir);
end

for isnr=1:numel(snrList)
    snrDb=snrList(isnr);sigma=sigmaRef(isnr);
    for seed=1:numSeeds
        done=true;
        for ist=1:numel(strategyKeys)
            done=done&&is_completed(rows,ir,seed,snrDb,strategyKeys(ist));
        end
        if done,continue;end
        lowSeed=971000+seed;highSeed=972000+seed;
        lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
            cfg0,lib.domain,sigma,lowSeed,'fixed_std');
        calTic=tic;C=calibrate_inv_log_2_low_speed(lowData,lib,cfg0);
        calibrationTime=toc(calTic);
        high=simulate_highspeed_from_low_increment(truthModel,deltaTrue,cfg0,...
            sigma,'fixed_std',highSeed);
        highMap=map_highspeed_to_space(high,cfg0.alpha_k,cfg0.R_tip,lib.domain,.02);
        for ist=1:numel(strategyKeys)
            if is_completed(rows,ir,seed,snrDb,strategyKeys(ist)),continue;end
            cfg=cfg0;cfg.dualSyncSearchStrategy=strategyKeys(ist);
            fitTic=tic;
            fit=run_inv_log_2_high_with_calibration(C,highMap,cfg);
            fitTime=toc(fitTic);
            ir=ir+1;
            rows(ir)=pack_row(seed,snrDb,sigma,lowSeed,highSeed,strategyKeys(ist),...
                strategyLabels(ist),C,fit,calibrationTime,fitTime,truthT,lib.xGrid,...
                gLow,gHigh,deltaTrue,cfg0.A_true,eoTrue,high);
            Detail=struct2table(rows(1:ir));writetable(Detail,checkpoint);
            save(fullfile(outDir,'checkpoint.mat'),'Detail','ReferenceNoise','cfg0',...
                'eoTrue','gLow','gHigh','deltaTrue','-v7.3');
            fprintf('SNR=%g dB seed=%d %s: EO=(%d,%d), max |dA|=%.4f mm, |ddg|=%.4g mm, %.1f s\n',...
                snrDb,seed,strategyLabels(ist),fit.eo_id(1),fit.eo_id(2),...
                rows(ir).max_amplitude_error_mm,rows(ir).delta_gap_error_mm,fitTime);
        end
    end
end

Detail=struct2table(rows(1:ir));Summary=summarize_detail(Detail,snrList,strategyKeys,strategyLabels);
writetable(Detail,fullfile(outDir,'detail.csv'));
writetable(Summary,fullfile(outDir,'summary.csv'));
save(fullfile(outDir,'paper_baseline_fixed_noise_eo_strategy.mat'),'Detail','Summary',...
    'ReferenceNoise','cfg0','eoTrue','gLow','gHigh','deltaTrue','-v7.3');
plot_results(Summary,snrList,strategyKeys,strategyLabels,outDir);
disp(Summary);
end

function sigma=reference_noise_std(truthModel,delta,cfg,domain,snrList)
clean=simulate_highspeed_from_low_increment(truthModel,delta,cfg,Inf,'snr_db',1);
cfgStatic=cfg;cfgStatic.A_true=[];cfgStatic.f_true=[];cfgStatic.phi_true=[];
stat=simulate_highspeed_from_low_increment(truthModel,delta,cfgStatic,Inf,'snr_db',1);
mv=map_highspeed_to_space(clean,cfg.alpha_k,cfg.R_tip,domain,.02);
ms=map_highspeed_to_space(stat,cfg.alpha_k,cfg.R_tip,domain,.02);
vibrationRms=rms(mv.V_a-ms.V_a);
sigma=vibrationRms./10.^(snrList/20);
end

function tf=is_completed(rows,n,seed,snrDb,key)
tf=false;
for i=1:n
    if rows(i).seed==seed&&rows(i).reference_snr_db==snrDb&&...
            string(rows(i).strategy_key)==key,tf=true;return;end
end
end

function r=empty_row()
r=struct('seed',NaN,'reference_snr_db',NaN,'fixed_noise_std_V',NaN,...
    'low_seed',NaN,'high_seed',NaN,'strategy_key',"",'strategy_label',"",...
    'selected_window_mm',NaN,'selected_span_bins',NaN,'template_rmse_V',NaN,...
    'weighted_derivative_rmse_V_per_mm',NaN,'g_low_true_mm',NaN,'g_low_est_mm',NaN,...
    'g0_error_mm',NaN,'g_high_true_mm',NaN,'g_high_est_mm',NaN,...
    'high_gap_error_mm',NaN,'delta_gap_true_mm',NaN,'delta_gap_est_mm',NaN,...
    'delta_gap_error_mm',NaN,'eo1_est',NaN,'eo2_est',NaN,'eo_correct',false,...
    'A1_est_mm',NaN,'A2_est_mm',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,...
    'max_amplitude_error_mm',NaN,'parameter_success',false,'rmse_V',NaN,...
    'true_pair_rank',NaN,'true_pair_refined',false,'true_pair_refine_start_count',NaN,...
    'calibration_time_s',NaN,'fit_time_s',NaN,'nonlinear_refine_count',NaN,...
    'low_equivalent_snr_db',NaN,'high_equivalent_snr_db',NaN);
end

function r=pack_row(seed,snrDb,sigma,lowSeed,highSeed,key,label,C,fit,calTime,fitTime,...
        truthT,x,gLow,gHigh,deltaTrue,Atrue,eoTrue,high)
r=empty_row();r.seed=seed;r.reference_snr_db=snrDb;r.fixed_noise_std_V=sigma;
r.low_seed=lowSeed;r.high_seed=highSeed;r.strategy_key=key;r.strategy_label=label;
r.selected_window_mm=get_field(C,'low_speed_selected_window_mm',NaN);
r.selected_span_bins=get_field(C,'low_speed_selected_span_bins',NaN);
T=C.lowTemplate.templateBySensor;d=median(diff(x));D=spatial_gradient(T,d);
Dt=spatial_gradient(truthT,d);w=Dt.^2;
r.template_rmse_V=sqrt(mean((T-truthT).^2,'all'));
r.weighted_derivative_rmse_V_per_mm=sqrt(sum(w.*(D-Dt).^2,'all')/sum(w,'all'));
r.g_low_true_mm=gLow;r.g_low_est_mm=C.pathCal.g0;r.g0_error_mm=abs(C.pathCal.g0-gLow);
r.g_high_true_mm=gHigh;r.g_high_est_mm=fit.g_used;r.high_gap_error_mm=abs(fit.g_used-gHigh);
r.delta_gap_true_mm=deltaTrue;r.delta_gap_est_mm=fit.g_used-C.pathCal.g0;
r.delta_gap_error_mm=abs(r.delta_gap_est_mm-deltaTrue);
eo=fit.eo_id(:).';A=fit.A_id(:).';r.eo1_est=eo(1);r.eo2_est=eo(2);
r.eo_correct=isequal(eo,eoTrue);r.A1_est_mm=A(1);r.A2_est_mm=A(2);
r.A1_error_mm=abs(A(1)-Atrue(1));r.A2_error_mm=abs(A(2)-Atrue(2));
r.max_amplitude_error_mm=max(r.A1_error_mm,r.A2_error_mm);
r.parameter_success=r.eo_correct&&r.max_amplitude_error_mm<=.01;
r.rmse_V=fit.rmse;r.calibration_time_s=calTime;r.fit_time_s=fitTime;
r.nonlinear_refine_count=get_field(fit,'nonlinear_refine_count',NaN);
r.low_equivalent_snr_db=NaN;r.high_equivalent_snr_db=high.snr_db_equiv;
if isfield(fit,'diagnostic_pair_eo')
    q=find(all(fit.diagnostic_pair_eo==eoTrue,2),1);
    if ~isempty(q),r.true_pair_rank=fit.diagnostic_pair_rank(q);end
end
if isfield(fit,'diagnostic_refine_start_eo')
    r.true_pair_refine_start_count=sum(all(fit.diagnostic_refine_start_eo==eoTrue,2));
end
if isfield(fit,'diagnostic_refined_eo')
    r.true_pair_refined=any(all(fit.diagnostic_refined_eo==eoTrue,2));
end
end

function S=summarize_detail(D,snrList,keys,labels)
s=repmat(empty_summary(),numel(snrList)*numel(keys),1);ir=0;
for isnr=1:numel(snrList)
    for ik=1:numel(keys)
        q=D.reference_snr_db==snrList(isnr)&D.strategy_key==keys(ik);X=D(q,:);ir=ir+1;
        s(ir).reference_snr_db=snrList(isnr);s(ir).strategy_key=keys(ik);
        s(ir).strategy_label=labels(ik);s(ir).n=height(X);
        s(ir).eo_correct_rate=mean(X.eo_correct);s(ir).parameter_success_rate=mean(X.parameter_success);
        s(ir).max_amplitude_error_mm_median=median(X.max_amplitude_error_mm);
        s(ir).max_amplitude_error_mm_max=max(X.max_amplitude_error_mm);
        s(ir).delta_gap_error_mm_median=median(X.delta_gap_error_mm);
        s(ir).delta_gap_error_mm_max=max(X.delta_gap_error_mm);
        s(ir).g0_error_mm_median=median(X.g0_error_mm);
        s(ir).selected_window_mm_median=median(X.selected_window_mm);
        s(ir).template_rmse_V_median=median(X.template_rmse_V);
        s(ir).derivative_rmse_median=median(X.weighted_derivative_rmse_V_per_mm);
        s(ir).fit_time_s_median=median(X.fit_time_s);
        s(ir).nonlinear_refine_count_median=median(X.nonlinear_refine_count,'omitnan');
    end
end
S=struct2table(s);
end

function s=empty_summary()
s=struct('reference_snr_db',NaN,'strategy_key',"",'strategy_label',"",'n',0,...
    'eo_correct_rate',NaN,'parameter_success_rate',NaN,...
    'max_amplitude_error_mm_median',NaN,'max_amplitude_error_mm_max',NaN,...
    'delta_gap_error_mm_median',NaN,'delta_gap_error_mm_max',NaN,...
    'g0_error_mm_median',NaN,'selected_window_mm_median',NaN,...
    'template_rmse_V_median',NaN,'derivative_rmse_median',NaN,...
    'fit_time_s_median',NaN,'nonlinear_refine_count_median',NaN);
end

function plot_results(S,snrList,keys,labels,outDir)
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 6.8]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');cc=[.30 .48 .70;.85 .45 .30];
metrics={'eo_correct_rate','max_amplitude_error_mm_median','fit_time_s_median'};
ylabels={'EO correct rate (%)','Median max amplitude error (\mum)','Median fit time (s)'};
scales=[100 1000 1];bh=[];
for ip=1:3
    ax=nexttile;Y=nan(numel(snrList),numel(keys));
    for i=1:numel(snrList)
        for k=1:numel(keys)
            q=S.reference_snr_db==snrList(i)&S.strategy_key==keys(k);
            Y(i,k)=S.(metrics{ip})(q)*scales(ip);
        end
    end
    b=bar(1:numel(snrList),Y,'grouped');for k=1:numel(b),b(k).FaceColor=cc(k,:);end
    if ip==1,bh=b;ylim([0 105]);end
    set(ax,'XTick',1:numel(snrList),'XTickLabel',string(snrList));xlabel('Reference SNR (dB)');ylabel(ylabels{ip});
    text(ax,-.13,1.04,char('a'+ip-1),'Units','normalized','FontWeight','bold');
end
set(findall(fig,'Type','axes'),'TickDir','in','Box','on','FontName','Times New Roman','FontSize',8,'LineWidth',.7);
lg=legend(bh,cellstr(labels),'Box','off','Orientation','horizontal');lg.Layout.Tile='south';
exportgraphics(fig,fullfile(outDir,'fixed_noise_eo_strategy_comparison.png'),'Resolution',300);
try
    exportgraphics(fig,fullfile(outDir,'fixed_noise_eo_strategy_comparison.pdf'),'ContentType','vector');
catch
end
close(fig);
end

function D=spatial_gradient(T,d)
if isvector(T),D=gradient(T,d);else,[~,D]=gradient(T,d,d);end
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
