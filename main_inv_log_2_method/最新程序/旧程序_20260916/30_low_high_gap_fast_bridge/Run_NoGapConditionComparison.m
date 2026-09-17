function T=Run_NoGapConditionComparison(seedOffset)
% Compare the gap-aware structured route on the four no-gap benchmark cases.
if nargin<1,seedOffset=0;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=configure_structured_identification(base.cfgCase,"experiment_current");
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=4;
cfg.structuredMinRevolutions=4;cfg.structuredRequireMinRevolutions=true;
cfg.route30ForwardModel="low_increment";cfg.route30LowTemplateMethod="adaptive_sg";
cfg.route30SupportAware=true;cfg.route30DualFrequencyRangeHz=[300 1500];
cfg.route30SingleFrequencyRangeHz=[300 1500];cfg.route30StructuredFastBudgets=false;
cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg.dualSyncGapCount=11;cfg.dualSyncKeepPerGap=80;cfg.dualSyncIteratedReplayCount=2100;cfg.dualSyncRefineCount=50;
cfg.syncAsyncGapCount=15;cfg.syncAsyncKeepPerGap=100;cfg.syncAsyncIteratedReplayCount=2100;cfg.syncAsyncRefineCount=50;
cfg.dualSyncSearchStrategy="funnel_v1";
cfg.funnelAllPairsFirstJointStep=true;cfg.funnelSecondJointPairCount=24;cfg.funnelFullRefineCount=3;
cfg.structuredComponentAmplitudeFloorMm=.05;
cases=[case_row("500+1300","dual_sync_sync",[500 1300],[10 26]);
       case_row("500+560","dual_sync_async",[500 560],[10 NaN]);
       case_row("700+1200","dual_sync_sync",[700 1200],[14 24]);
       case_row("900+1400","dual_sync_sync",[900 1400],[18 28])];
snrs=[5 10 15 20];rows=repmat(row0(),numel(cases)*numel(snrs),1);ir=0;
for isnr=1:numel(snrs)
  snr=snrs(isnr);
  for ic=1:numel(cases)
    c=cases(ic);cfg.snrDb=snr;cfg.route30FrequencyStructureMode=c.mode;
    cfg.f_true=c.f;cfg.A_true=[.25 .15];cfg.phi_true=[pi/4 -pi/6];
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,.5,z),cfg,lib.domain,Inf,710000+100*ic+snr+seedOffset,'fixed_std');
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,.6,z),cfg.RPM_high,4,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,[.25 .15],c.f,[pi/4 -pi/6],'fixed_std');
    sigma=sqrt(var(d.V_clean,1))/10^(snr/20);rng(730000+100*ic+snr+seedOffset,'twister');d.V_cap=d.V_clean+sigma*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    ticFit=tic;fit=run_inv_log_2_structured_main(low,map,lib,cfg,c.mode);elapsed=toc(ticFit);
    fEst=sort(fit.f_id(:).');fErr=max(abs(fEst-sort(c.f)));
    if c.mode=="dual_sync_async",eoOk=fit.sync_eo==10&&abs(fit.async_frequency_hz-560)<=1;else,eoOk=isequal(fit.eo_id,c.eo);end
    ir=ir+1;rows(ir)=struct('case_id',c.id,'mode',c.mode,'snr_db',snr,'f_true',string(mat2str(c.f,7)),...
      'f_est',string(mat2str(fEst,7)),'max_frequency_error_hz',fErr,'eo_or_structure_correct',eoOk,...
      'amplitude_est',string(mat2str(fit.A_id,7)),'g_true',.6,'g_est',fit.g_used,...
      'gap_error_mm',abs(fit.g_used-.6),'rmse_v',fit.rmse,'identification_status',fit.identification_status,...
      'identification_confident',fit.identification_confident,'elapsed_s',elapsed,'success',eoOk&&fErr<=1&&abs(fit.g_used-.6)<=.05);
  end
end
T=struct2table(rows);T.seed_offset=repmat(seedOffset,height(T),1);out=fullfile(root,'output','no_gap_condition_comparison');if ~exist(out,'dir'),mkdir(out);end
writetable(T,fullfile(out,sprintf('gap_aware_seed_%d.csv',seedOffset)));save(fullfile(out,sprintf('gap_aware_seed_%d.mat',seedOffset)),'T');disp(T);
end
function c=case_row(id,mode,f,eo),c=struct('id',id,'mode',mode,'f',f,'eo',eo);end
function r=row0(),r=struct('case_id',"",'mode',"",'snr_db',NaN,'f_true',"",'f_est',"",'max_frequency_error_hz',NaN,'eo_or_structure_correct',false,'amplitude_est',"",'g_true',NaN,'g_est',NaN,'gap_error_mm',NaN,'rmse_v',NaN,'identification_status',"",'identification_confident',false,'elapsed_s',NaN,'success',false);end
