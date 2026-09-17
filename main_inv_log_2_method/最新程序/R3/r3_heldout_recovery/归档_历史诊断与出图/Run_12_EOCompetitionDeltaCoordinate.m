function Result = Run_12_EOCompetitionDeltaCoordinate()
%RUN_12_EOCOMPETITIONDELTACOORDINATE Recheck EO failures at consistent state.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
S=load(fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat'),'Result');
C0=S.Result.C0;sigmaV=S.Result.sigmaV;
Formal=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
Fail=Formal(Formal.method=="adaptive"&Formal.A_true_mm>0&Formal.vib_success==0,:);
assert(height(Fail)==8,'Expected the eight formal adaptive nonzero failures.');
rows=repmat(empty_row(),height(Fail),1);
for i=1:height(Fail)
    g=Fail.g_truth_mm(i);A=Fail.A_true_mm(i);seed=Fail.high_seed(i);
    rep=find(P.mainHighSeeds==seed,1);
    [~,highMap]=r3_generate_high_case(R3,g,A,P.phaseValuesRad(rep),sigmaV,seed);
    qDelta=C0.pathCal.g0+(g-P.gReferenceMm);
    joint10=fit_eo(R3,C0,highMap,P,10,false,qDelta);
    joint16=fit_eo(R3,C0,highMap,P,16,false,qDelta);
    delta10=fit_eo(R3,C0,highMap,P,10,true,qDelta);
    delta16=fit_eo(R3,C0,highMap,P,16,true,qDelta);
    rows(i)=pack(Fail(i,:),rep,qDelta,joint10,joint16,delta10,delta16);
    fprintf('Delta EO diagnostic %d/%d: g=%.3f A=%.3f seed=%d\n',i,height(Fail),g,A,seed);
end
Detail=struct2table(rows);
outDir=fullfile(root,'output','12_eo_competition_delta_coordinate');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Detail,fullfile(outDir,'eo_competition_delta_detail.csv'));
Summary=groupcounts(Detail,'classification');
writetable(Summary,fullfile(outDir,'eo_competition_delta_summary.csv'));
Result=struct('P',P,'Detail',Detail,'Summary',Summary,'outputDir',outDir);
save(fullfile(outDir,'eo_competition_delta_result.mat'),'Result','-v7.3');
fprintf('EO_DELTA_DIAGNOSTIC_OK rows=%d\n',height(Detail));
end

function fit=fit_eo(R3,C0,highMap,P,eo,fixed,qDelta)
cfg=R3.cfg;cfg.singleSyncCandidateEO=eo;cfg.syncSingleKeepPerGap=1;cfg.syncSingleRefineCount=1;
state=C0.lowState;
if fixed
    cfg.route30GapHalfWidthMm=0;cfg.route30GapCount=1;cfg.route30SingleGapHalfWidthMm=0;state.gHat=qDelta;
else
    cfg.route30GapHalfWidthMm=0.5*diff(P.adaptiveSearchBoundsMm);
    cfg.route30GapCount=numel(P.adaptiveCoarseGapGridMm);cfg.route30GapSearchBoundsMm=P.adaptiveSearchBoundsMm;
    cfg.route30SingleGapHalfWidthMm=.02;state.gHat=mean(P.adaptiveSearchBoundsMm);
end
fit=run_single_sync_voltage_vp(highMap,C0.templateModel,cfg,state);
end

function r=pack(F,rep,q,j10,j16,d10,d16)
r=empty_row();r.g_truth_mm=F.g_truth_mm;r.A_true_mm=F.A_true_mm;r.high_seed=F.high_seed;
r.replicate_index=rep;r.formal_adaptive_eo=F.eo_est;r.formal_adaptive_g_mm=F.g_est_mm;
r.delta_coordinate_mm=q;r.joint_eo10_rmse_V=j10.rmse;r.joint_eo16_rmse_V=j16.rmse;
r.delta_eo10_rmse_V=d10.rmse;r.delta_eo16_rmse_V=d16.rmse;
r.joint_delta_16_minus_10_V=j16.rmse-j10.rmse;
r.fixed_delta_16_minus_10_V=d16.rmse-d10.rmse;
if r.joint_delta_16_minus_10_V<0 && r.fixed_delta_16_minus_10_V>0
    r.classification="joint_search_rank_inversion_confirmed";
elseif r.joint_delta_16_minus_10_V<0 && r.fixed_delta_16_minus_10_V<0
    r.classification="eo16_better_even_at_delta_state";
else
    r.classification="joint_eo16_preference_not_reproduced";
end
end

function r=empty_row()
r=struct('g_truth_mm',NaN,'A_true_mm',NaN,'high_seed',NaN,'replicate_index',NaN,...
    'formal_adaptive_eo',NaN,'formal_adaptive_g_mm',NaN,'delta_coordinate_mm',NaN,...
    'joint_eo10_rmse_V',NaN,'joint_eo16_rmse_V',NaN,'delta_eo10_rmse_V',NaN,...
    'delta_eo16_rmse_V',NaN,'joint_delta_16_minus_10_V',NaN,...
    'fixed_delta_16_minus_10_V',NaN,'classification',"");
end
