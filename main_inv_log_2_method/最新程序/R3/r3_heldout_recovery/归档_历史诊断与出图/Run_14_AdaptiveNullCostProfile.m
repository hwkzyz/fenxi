function Result = Run_14_AdaptiveNullCostProfile()
%RUN_14_ADAPTIVENULLCOSTPROFILE Profile null and vibration cost versus gap.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
S=load(fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat'),'Result');
C0=S.Result.C0;sigmaV=S.Result.sigmaV;
T=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');
T=T(T.method=="adaptive"&T.A_true_mm==0,:);

selected=table();
for g=[0.5 0.8 0.9 1.0 1.2]
    Q=T(T.g_truth_mm==g&T.false_positive==1,:);assert(~isempty(Q));selected=[selected;Q(1,:)]; %#ok<AGROW>
end
Q=T(T.g_truth_mm==0.8&T.false_positive==0,:);assert(~isempty(Q));selected=[selected;Q(1,:)];

rows=repmat(empty_row(),height(selected)*21,1);k=0;
for i=1:height(selected)
    F=selected(i,:);g=F.g_truth_mm;seed=F.high_seed;
    rep=find(P.mainHighSeeds==seed,1);[~,map]=r3_generate_high_case(R3,g,0,P.phaseValuesRad(rep),sigmaV,seed);
    qExpected=C0.pathCal.g0+(g-P.gReferenceMm);
    qGrid=unique([linspace(P.adaptiveSearchBoundsMm(1),P.adaptiveSearchBoundsMm(2),19),...
        qExpected,F.g_model_coordinate_mm]);
    for q=qGrid
        fit=fit_fixed_gap(R3,C0,map,q);
        [nullRmse,nullDx]=fit_null_dx(map,C0.templateModel,R3.cfg,q,C0.lowState.dx0);
        k=k+1;rows(k)=pack(i,F,q,qExpected,fit,nullRmse,nullDx);
    end
    fprintf('Null profile %d/%d: g=%.3f seed=%d FP=%d\n',i,height(selected),g,seed,F.false_positive);
end
Detail=struct2table(rows(1:k));Summary=repmat(empty_summary(),height(selected),1);
for i=1:height(selected)
    Q=Detail(Detail.case_index==i,:);[~,j]=min(Q.vibration_rmse_V);B=Q(j,:);
    [~,j0]=min(Q.null_rmse_V);N=Q(j0,:);
    E=Q(abs(Q.g_query_mm-Q.expected_coordinate_mm)<1e-10,:);assert(height(E)==1);
    Summary(i)=pack_summary(B,N,E);
end
Summary=struct2table(Summary);
outDir=fullfile(root,'output','14_adaptive_null_cost_profile');if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Detail,fullfile(outDir,'adaptive_null_cost_profile_detail.csv'));
writetable(Summary,fullfile(outDir,'adaptive_null_cost_profile_summary.csv'));
Result=struct('P',P,'Detail',Detail,'Summary',Summary,'outputDir',outDir);
save(fullfile(outDir,'adaptive_null_cost_profile.mat'),'Result','-v7.3');
fprintf('ADAPTIVE_NULL_COST_PROFILE_OK cases=%d rows=%d\n',height(Summary),height(Detail));
end

function fit=fit_fixed_gap(R3,C0,map,q)
cfg=R3.cfg;cfg.route30GapHalfWidthMm=0;cfg.route30GapCount=1;cfg.route30SingleGapHalfWidthMm=0;
state=C0.lowState;state.gHat=q;fit=run_single_sync_voltage_vp(map,C0.templateModel,cfg,state);
end

function [rmse,dx]=fit_null_dx(map,lib,cfg,q,dx0)
half=.5;if isfield(cfg,'route30DxHalfWidthMm'),half=cfg.route30DxHalfWidthMm;end
obj=@(z)null_mse(z,map,lib,q);opts=optimset('Display','off','TolX',1e-9);
dx=fminbnd(obj,dx0-half,dx0+half,opts);rmse=sqrt(obj(dx));
end

function mse=null_mse(dx,map,lib,q)
r=map.V_a(:)-eval_gap_template(lib,q,map.x_v(:)-dx);r=r(isfinite(r));mse=mean(r.^2);
end

function r=pack(i,F,q,qExpected,fit,nullRmse,nullDx)
r=empty_row();r.case_index=i;r.g_truth_mm=F.g_truth_mm;r.high_seed=F.high_seed;
r.formal_false_positive=F.false_positive;r.formal_g_query_mm=F.g_model_coordinate_mm;
r.formal_A_est_mm=F.A_est_mm;r.g_query_mm=q;r.expected_coordinate_mm=qExpected;
r.coordinate_error_mm=q-qExpected;r.null_rmse_V=nullRmse;r.null_dx_mm=nullDx;
r.vibration_rmse_V=fit.rmse;r.rmse_improvement_V=nullRmse-fit.rmse;
r.eo_est=fit.eo_id;r.A_est_mm=fit.A_id;r.exceeds_threshold=fit.A_id>F.detection_threshold_mm;
end

function s=pack_summary(B,N,E)
s=empty_summary();s.case_index=B.case_index;s.g_truth_mm=B.g_truth_mm;s.high_seed=B.high_seed;
s.formal_false_positive=B.formal_false_positive;s.expected_coordinate_mm=B.expected_coordinate_mm;
s.best_vibration_g_mm=B.g_query_mm;s.best_vibration_coordinate_error_mm=B.coordinate_error_mm;
s.best_vibration_eo=B.eo_est;s.best_vibration_A_mm=B.A_est_mm;
s.best_vibration_rmse_V=B.vibration_rmse_V;s.best_null_g_mm=N.g_query_mm;s.best_null_rmse_V=N.null_rmse_V;
s.expected_state_A_mm=E.A_est_mm;s.expected_state_eo=E.eo_est;
s.expected_state_rmse_improvement_V=E.rmse_improvement_V;s.expected_state_exceeds_threshold=E.exceeds_threshold;
end

function r=empty_row()
r=struct('case_index',NaN,'g_truth_mm',NaN,'high_seed',NaN,'formal_false_positive',false,...
    'formal_g_query_mm',NaN,'formal_A_est_mm',NaN,'g_query_mm',NaN,'expected_coordinate_mm',NaN,...
    'coordinate_error_mm',NaN,'null_rmse_V',NaN,'null_dx_mm',NaN,'vibration_rmse_V',NaN,...
    'rmse_improvement_V',NaN,'eo_est',NaN,'A_est_mm',NaN,'exceeds_threshold',false);
end

function s=empty_summary()
s=struct('case_index',NaN,'g_truth_mm',NaN,'high_seed',NaN,'formal_false_positive',false,...
    'expected_coordinate_mm',NaN,'best_vibration_g_mm',NaN,'best_vibration_coordinate_error_mm',NaN,...
    'best_vibration_eo',NaN,'best_vibration_A_mm',NaN,'best_vibration_rmse_V',NaN,...
    'best_null_g_mm',NaN,'best_null_rmse_V',NaN,'expected_state_A_mm',NaN,'expected_state_eo',NaN,...
    'expected_state_rmse_improvement_V',NaN,'expected_state_exceeds_threshold',false);
end
