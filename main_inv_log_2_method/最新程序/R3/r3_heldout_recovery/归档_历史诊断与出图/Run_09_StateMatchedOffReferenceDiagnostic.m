function Result = Run_09_StateMatchedOffReferenceDiagnostic(g)
%RUN_09_STATEMATCHEDOFFREFERENCEDIAGNOSTIC Compare state-match coordinates.
% Replays one formal gap using the formal records and writes an isolated output.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
if nargin<1||isempty(g),error('r3:MissingGap','Provide one formal R3 gap.');end
if ~ismembertol(g,P.gMainMm,1e-12,'DataScale',1),error('r3:UnknownGap','Gap is outside the formal matrix.');end
calFile=fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat');
detailFile=fullfile(root,'output','05_main_15db','main_detail_merged.csv');
S=load(calFile,'Result');C0=S.Result.C0;sigmaV=S.Result.sigmaV;
Formal=readtable(detailFile,'TextType','string');Aall=[P.amplitudeNullMm P.amplitudeMainMm];
rows=repmat(empty_row(),numel(Aall)*numel(P.mainHighSeeds)*2,1);k=0;
for A=Aall
    for rep=1:numel(P.mainHighSeeds)
        seed=P.mainHighSeeds(rep);phi=P.phaseValuesRad(rep);
        [~,highMap,truth]=r3_generate_high_case(R3,g,A,phi,sigmaV,seed);
        formal=Formal(Formal.method=="state_matched"&Formal.g_truth_mm==g&...
            Formal.A_true_mm==A&Formal.high_seed==seed,:);
        if height(formal)~=1,error('r3:FormalRow','Expected one formal SM row.');end
        qAbs=g;qDelta=C0.pathCal.g0+(g-P.gReferenceMm);
        k=k+1;rows(k)=pack("SM_abs",qAbs,fit_one(R3,C0,highMap,qAbs),truth,P,C0,formal,rep);
        k=k+1;rows(k)=pack("SM_delta",qDelta,fit_one(R3,C0,highMap,qDelta),truth,P,C0,formal,rep);
    end
    fprintf('Completed g=%.3f mm, A=%.3f mm\n',g,A);
end
Detail=struct2table(rows(1:k));
Summary=groupsummary(Detail,{'coordinate_mode','A_true_mm'},'mean',...
    {'vib_success','false_positive','reproduces_formal_sm'});
Summary.Properties.VariableNames{'mean_vib_success'}='P_vib';
Summary.Properties.VariableNames{'mean_false_positive'}='P_FP';
Summary.Properties.VariableNames{'mean_reproduces_formal_sm'}='P_reproduces_formal_sm';
outDir=fullfile(root,'output','09_state_matched_offreference_diagnostic');if ~exist(outDir,'dir'),mkdir(outDir);end
tag=sprintf('g%04d',round(1000*g));
writetable(Detail,fullfile(outDir,[tag '_detail.csv']));writetable(Summary,fullfile(outDir,[tag '_summary.csv']));
Result=struct('P',P,'C0',C0,'Detail',Detail,'Summary',Summary,'outputDir',outDir);
save(fullfile(outDir,[tag '_result.mat']),'Result','-v7.3');
if any(Detail.coordinate_mode=="SM_abs"&~Detail.reproduces_formal_sm),error('r3:NonreproducibleSM','SM_abs replay mismatch.');end
fprintf('SM_OFFREFERENCE_OK g=%.3f rows=%d\n',g,height(Detail));
end

function fit=fit_one(R3,C0,map,g)
cfg=R3.cfg;cfg.route30GapHalfWidthMm=0;cfg.route30GapCount=1;cfg.route30SingleGapHalfWidthMm=0;
state=C0.lowState;state.gHat=g;fit=run_single_sync_voltage_vp(map,C0.templateModel,cfg,state);
end

function r=pack(mode,gq,fit,truth,P,C0,formal,rep)
ampErr=fit.A_id-truth.A_truth_mapped_mm;phaseErr=NaN;success=false;
if truth.A_true_mm>0
    phaseErr=abs(angle(exp(1i*(fit.phi_id-truth.phi_mapped_rad))));
    tol=max(P.amplitudeAbsToleranceMm,P.amplitudeRelTolerance*truth.A_truth_mapped_mm);
    success=fit.eo_id==P.eoTrue&&abs(ampErr)<=tol&&phaseErr<=P.phaseToleranceRad;
end
r=empty_row();r.coordinate_mode=mode;r.g_truth_mm=truth.g_truth_mm;r.A_true_mm=truth.A_true_mm;
r.replicate_index=rep;r.high_seed=truth.high_seed;r.g_query_mm=gq;
r.increment_from_g0_mm=gq-C0.pathCal.g0;r.true_increment_mm=truth.g_truth_mm-P.gReferenceMm;
r.increment_error_mm=r.increment_from_g0_mm-r.true_increment_mm;r.eo_est=fit.eo_id;r.A_est_mm=fit.A_id;
r.phase_error_rad=phaseErr;r.rmse_V=fit.rmse;r.vib_success=success;
r.false_positive=truth.A_true_mm==0&&fit.A_id>formal.detection_threshold_mm;
r.reproduces_formal_sm=fit.eo_id==formal.eo_est&&abs(fit.A_id-formal.A_est_mm)<1e-10&&abs(fit.rmse-formal.rmse_V)<1e-12;
end

function r=empty_row()
r=struct('coordinate_mode',"",'g_truth_mm',NaN,'A_true_mm',NaN,'replicate_index',NaN,...
    'high_seed',NaN,'g_query_mm',NaN,'increment_from_g0_mm',NaN,'true_increment_mm',NaN,...
    'increment_error_mm',NaN,'eo_est',NaN,'A_est_mm',NaN,'phase_error_rad',NaN,'rmse_V',NaN,...
    'vib_success',false,'false_positive',false,'reproduces_formal_sm',false);
end
