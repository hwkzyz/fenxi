function Result = Run_06_EOCompetitionDiagnostic()
%RUN_06_EOCOMPETITIONDIAGNOSTIC Replay failed adaptive fits by fixed EO.
% This diagnostic only reads formal R3 results and writes a separate output.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
calFile=fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat');
detailFile=fullfile(root,'output','05_main_15db','main_detail_merged.csv');
if ~exist(calFile,'file')||~exist(detailFile,'file'),error('r3:MissingFormalInput','Missing frozen calibration or merged formal R3 detail.');end
S=load(calFile,'Result');C0=S.Result.C0;sigmaV=S.Result.sigmaV;
Detail=readtable(detailFile,'TextType','string');
q=Detail.method=="adaptive"&Detail.A_true_mm>0&Detail.vib_success==0;
Fail=Detail(q,:);
if isempty(Fail),error('r3:NoAdaptiveFailures','No adaptive non-null failures were found.');end

outDir=fullfile(root,'output','06_eo_competition_diagnostic');
if ~exist(outDir,'dir'),mkdir(outDir);end
rows=repmat(empty_row(),height(Fail),1);
for i=1:height(Fail)
    g=Fail.g_truth_mm(i);A=Fail.A_true_mm(i);seed=Fail.high_seed(i);
    rep=find(P.mainHighSeeds==seed,1);
    if isempty(rep),error('r3:UnknownSeed','Seed %d is not in P.mainHighSeeds.',seed);end
    [~,highMap]=r3_generate_high_case(R3,g,A,P.phaseValuesRad(rep),sigmaV,seed);
    joint10=fit_one_eo(R3,C0,highMap,P,10,false,g);
    joint16=fit_one_eo(R3,C0,highMap,P,16,false,g);
    matched10=fit_one_eo(R3,C0,highMap,P,10,true,g);
    matched16=fit_one_eo(R3,C0,highMap,P,16,true,g);
    rows(i)=pack_row(g,A,seed,rep,P.phaseValuesRad(rep),Fail.eo_est(i),Fail.g_est_mm(i),...
        joint10,joint16,matched10,matched16);
    fprintf('Diagnostic %d/%d: g=%.3f, A=%.3f, seed=%d\n',i,height(Fail),g,A,seed);
end
Diagnostic=struct2table(rows);
writetable(Diagnostic,fullfile(outDir,'eo_competition_detail.csv'));
Result=struct('P',P,'Diagnostic',Diagnostic,'sourceDetailFile',detailFile,'outputDir',outDir);
save(fullfile(outDir,'eo_competition_result.mat'),'Result','-v7.3');
fprintf('EO_DIAGNOSTIC_OK rows=%d output=%s\n',height(Diagnostic),outDir);
end

function fit=fit_one_eo(R3,C0,highMap,P,eo,fixGap,gTruth)
cfg=R3.cfg;cfg.singleSyncCandidateEO=eo;cfg.syncSingleKeepPerGap=1;cfg.syncSingleRefineCount=1;
state=C0.lowState;
if fixGap
    cfg.route30GapHalfWidthMm=0;cfg.route30GapCount=1;cfg.route30SingleGapHalfWidthMm=0;
    state.gHat=gTruth;
else
    cfg.route30GapHalfWidthMm=0.5*diff(P.adaptiveSearchBoundsMm);
    cfg.route30GapCount=numel(P.adaptiveCoarseGapGridMm);
    cfg.route30GapSearchBoundsMm=P.adaptiveSearchBoundsMm;
    cfg.route30SingleGapHalfWidthMm=.02;state.gHat=mean(P.adaptiveSearchBoundsMm);
end
fit=run_single_sync_voltage_vp(highMap,C0.templateModel,cfg,state);
end

function row=pack_row(g,A,seed,rep,phase,adaptiveEO,adaptiveG,j10,j16,m10,m16)
row=empty_row();row.g_truth_mm=g;row.A_true_mm=A;row.high_seed=seed;row.replicate_index=rep;
row.phase_injected_rad=phase;row.phase_injected_deg=phase*180/pi;
row.formal_adaptive_eo=adaptiveEO;row.formal_adaptive_g_mm=adaptiveG;
row.joint_eo10_rmse_V=j10.rmse;row.joint_eo10_g_mm=j10.g_used;row.joint_eo10_A_mm=j10.A_id;
row.joint_eo16_rmse_V=j16.rmse;row.joint_eo16_g_mm=j16.g_used;row.joint_eo16_A_mm=j16.A_id;
row.matched_eo10_rmse_V=m10.rmse;row.matched_eo16_rmse_V=m16.rmse;
row.delta_rmse_joint_16_minus_10_V=j16.rmse-j10.rmse;
row.delta_rmse_matched_16_minus_10_V=m16.rmse-m10.rmse;
if row.delta_rmse_joint_16_minus_10_V<0&&row.delta_rmse_matched_16_minus_10_V>0
    row.classification="joint_gap_rank_inversion";
elseif row.delta_rmse_joint_16_minus_10_V<0&&row.delta_rmse_matched_16_minus_10_V<0
    row.classification="eo16_competes_even_at_true_gap";
else
    row.classification="formal_failure_not_reproduced_by_fixed_eo_replay";
end
end

function row=empty_row()
row=struct('g_truth_mm',NaN,'A_true_mm',NaN,'high_seed',NaN,'replicate_index',NaN,...
    'phase_injected_rad',NaN,'phase_injected_deg',NaN,'formal_adaptive_eo',NaN,'formal_adaptive_g_mm',NaN,...
    'joint_eo10_rmse_V',NaN,'joint_eo10_g_mm',NaN,'joint_eo10_A_mm',NaN,...
    'joint_eo16_rmse_V',NaN,'joint_eo16_g_mm',NaN,'joint_eo16_A_mm',NaN,...
    'matched_eo10_rmse_V',NaN,'matched_eo16_rmse_V',NaN,...
    'delta_rmse_joint_16_minus_10_V',NaN,'delta_rmse_matched_16_minus_10_V',NaN,'classification',"");
end
