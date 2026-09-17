function Result = Run_08_StateMatchedCoordinateDiagnostic(amplitudes)
%RUN_08_STATEMATCHEDCOORDINATEDIAGNOSTIC Compare absolute and delta-matched queries.
% This diagnostic replays only the physical reference state (g=0.8 mm).

if nargin<1||isempty(amplitudes),amplitudes=R3_Protocol().amplitudeMainMm;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
if any(~ismember(amplitudes,[P.amplitudeNullMm P.amplitudeMainMm]))
    error('r3:UnknownAmplitude','Requested amplitude is outside the formal R3 grid.');
end
calFile=fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat');
detailFile=fullfile(root,'output','05_main_15db','main_detail_merged.csv');
if ~exist(calFile,'file')||~exist(detailFile,'file'),error('r3:MissingFormalInput','Missing frozen calibration or formal detail.');end
S=load(calFile,'Result');C0=S.Result.C0;sigmaV=S.Result.sigmaV;
Formal=readtable(detailFile,'TextType','string');
g=P.gReferenceMm;amplitudes=unique(amplitudes(:).','stable');
rows=repmat(empty_row(),numel(amplitudes)*numel(P.mainHighSeeds)*2,1);k=0;

for A=amplitudes
    for rep=1:numel(P.mainHighSeeds)
        seed=P.mainHighSeeds(rep);phi=P.phaseValuesRad(rep);
        [~,highMap,truth]=r3_generate_high_case(R3,g,A,phi,sigmaV,seed);
        absFit=fit_at_query(R3,C0,highMap,g);
        deltaQuery=C0.pathCal.g0+(g-P.gReferenceMm);
        deltaFit=fit_at_query(R3,C0,highMap,deltaQuery);
        formal=Formal(Formal.method=="state_matched"&Formal.g_truth_mm==g&...
            Formal.A_true_mm==A&Formal.high_seed==seed,:);
        if height(formal)~=1,error('r3:FormalRow','Expected exactly one formal state-matched row.');end
        k=k+1;rows(k)=pack_row("SM_abs",g,A,rep,seed,phi,g,absFit,truth,P,formal);
        k=k+1;rows(k)=pack_row("SM_delta",g,A,rep,seed,phi,deltaQuery,deltaFit,truth,P,formal);
    end
    fprintf('Completed g=%.3f mm, A=%.3f mm\n',g,A);
end
Detail=struct2table(rows(1:k));
Summary=groupsummary(Detail,{'coordinate_mode','A_true_mm'},'mean',...
    {'vib_success','joint_success','false_positive','reproduces_formal_sm'});
Summary.Properties.VariableNames{'mean_vib_success'}='P_vib';
Summary.Properties.VariableNames{'mean_joint_success'}='P_joint';
Summary.Properties.VariableNames{'mean_false_positive'}='P_FP';
Summary.Properties.VariableNames{'mean_reproduces_formal_sm'}='P_reproduces_formal_sm';
outDir=fullfile(root,'output','08_state_matched_coordinate_diagnostic');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Detail,fullfile(outDir,'reference_coordinate_detail.csv'));
writetable(Summary,fullfile(outDir,'reference_coordinate_summary.csv'));
Result=struct('P',P,'C0',C0,'Detail',Detail,'Summary',Summary,'outputDir',outDir);
save(fullfile(outDir,'reference_coordinate_result.mat'),'Result','-v7.3');
if any(Detail.coordinate_mode=="SM_abs"&~Detail.reproduces_formal_sm)
    error('r3:NonreproducibleSM','SM_abs replay did not reproduce formal state-matched results.');
end
fprintf('STATE_MATCHED_COORDINATE_OK rows=%d output=%s\n',height(Detail),outDir);
end

function fit=fit_at_query(R3,C0,highMap,gQuery)
cfg=R3.cfg;cfg.route30GapHalfWidthMm=0;cfg.route30GapCount=1;cfg.route30SingleGapHalfWidthMm=0;
state=C0.lowState;state.gHat=gQuery;
fit=run_single_sync_voltage_vp(highMap,C0.templateModel,cfg,state);
end

function row=pack_row(mode,g,A,rep,seed,phi,gQuery,fit,truth,P,formal)
ampErr=fit.A_id-truth.A_truth_mapped_mm;
if A==0
    phaseErr=NaN;vib=false;joint=false;
else
    phaseErr=abs(angle(exp(1i*(fit.phi_id-truth.phi_mapped_rad))));
    vib=fit.eo_id==P.eoTrue&&abs(ampErr)<=max(P.amplitudeAbsToleranceMm,P.amplitudeRelTolerance*A)&&phaseErr<=P.phaseToleranceRad;
    joint=vib&&abs(gQuery-g)<=P.gapToleranceMm;
end
row=empty_row();row.coordinate_mode=mode;row.g_truth_mm=g;row.A_true_mm=A;row.replicate_index=rep;
row.high_seed=seed;row.phase_injected_deg=phi*180/pi;row.g_query_mm=gQuery;
row.g_query_increment_from_g0_mm=gQuery-formal.g_low_model_mm;
row.eo_est=fit.eo_id;row.A_est_mm=fit.A_id;row.phase_error_rad=phaseErr;row.rmse_V=fit.rmse;
row.vib_success=vib;row.joint_success=joint;row.false_positive=(A==0&&fit.A_id>formal.detection_threshold_mm);
row.reproduces_formal_sm=fit.eo_id==formal.eo_est&&abs(fit.A_id-formal.A_est_mm)<1e-10&&abs(fit.rmse-formal.rmse_V)<1e-12;
end

function row=empty_row()
row=struct('coordinate_mode',"",'g_truth_mm',NaN,'A_true_mm',NaN,'replicate_index',NaN,...
    'high_seed',NaN,'phase_injected_deg',NaN,'g_query_mm',NaN,'g_query_increment_from_g0_mm',NaN,...
    'eo_est',NaN,'A_est_mm',NaN,'phase_error_rad',NaN,'rmse_V',NaN,'vib_success',false,...
    'joint_success',false,'false_positive',false,'reproduces_formal_sm',false);
end
