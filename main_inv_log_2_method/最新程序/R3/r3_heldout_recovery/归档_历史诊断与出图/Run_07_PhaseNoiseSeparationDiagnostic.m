function Result = Run_07_PhaseNoiseSeparationDiagnostic()
%RUN_07_PHASENOISESEPARATIONDIAGNOSTIC Separate phase and noise triggers.
% The formal R3 data couple phase index and seed. This replay changes one
% factor at a time for three representative adaptive EO=16 failures.

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
calFile=fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat');
if ~exist(calFile,'file'),error('r3:MissingC0','Missing frozen reference calibration.');end
S=load(calFile,'Result');C0=S.Result.C0;sigmaV=S.Result.sigmaV;

% One representative from each failed gap-amplitude configuration.
Cases=table([.9;.9;1.0],[.20;.25;.20],[2026082672;2026082674;2026082672],...
    [72;74;72],'VariableNames',{'g_truth_mm','A_true_mm','failure_seed','failure_rep'});
outDir=fullfile(root,'output','07_phase_noise_separation_diagnostic');
if ~exist(outDir,'dir'),mkdir(outDir);end
rows=repmat(empty_row(),2*75*height(Cases),1);k=0;
for ic=1:height(Cases)
    g=Cases.g_truth_mm(ic);A=Cases.A_true_mm(ic);seed0=Cases.failure_seed(ic);rep0=Cases.failure_rep(ic);
    for rep=1:75
        k=k+1;fit=fit_adaptive(R3,C0,g,A,P.phaseValuesRad(rep),sigmaV,seed0);
        rows(k)=pack_row(ic,"phase_sweep_fixed_noise",g,A,seed0,rep,P.phaseValuesRad(rep),fit);
    end
    for rep=1:75
        seed=P.mainHighSeeds(rep);fit=fit_adaptive(R3,C0,g,A,P.phaseValuesRad(rep0),sigmaV,seed);
        k=k+1;rows(k)=pack_row(ic,"seed_sweep_fixed_phase",g,A,seed,rep0,P.phaseValuesRad(rep0),fit);
    end
    fprintf('Completed representative case %d/%d: g=%.3f, A=%.3f\n',ic,height(Cases),g,A);
end
Detail=struct2table(rows);
Summary=groupsummary(Detail,{'case_id','axis','g_truth_mm','A_true_mm'},'mean','eo_success');
Summary.Properties.VariableNames{'mean_eo_success'}='P_EO10';
writetable(Detail,fullfile(outDir,'phase_noise_separation_detail.csv'));
writetable(Summary,fullfile(outDir,'phase_noise_separation_summary.csv'));
Result=struct('P',P,'Cases',Cases,'Detail',Detail,'Summary',Summary,'outputDir',outDir);
save(fullfile(outDir,'phase_noise_separation_result.mat'),'Result','-v7.3');
fprintf('PHASE_NOISE_DIAGNOSTIC_OK rows=%d output=%s\n',height(Detail),outDir);
end

function fit=fit_adaptive(R3,C0,g,A,phi,sigmaV,seed)
[~,highMap]=r3_generate_high_case(R3,g,A,phi,sigmaV,seed);
cfg=R3.cfg;P=R3.P;cfg.route30GapHalfWidthMm=.5*diff(P.adaptiveSearchBoundsMm);
cfg.route30GapCount=numel(P.adaptiveCoarseGapGridMm);cfg.route30GapSearchBoundsMm=P.adaptiveSearchBoundsMm;
cfg.route30SingleGapHalfWidthMm=.02;cfg.syncSingleKeepPerGap=25;
state=C0.lowState;state.gHat=mean(P.adaptiveSearchBoundsMm);
fit=run_single_sync_voltage_vp(highMap,C0.templateModel,cfg,state);
end

function row=pack_row(caseId,axis,g,A,seed,rep,phase,fit)
row=empty_row();row.case_id=caseId;row.axis=axis;row.g_truth_mm=g;row.A_true_mm=A;
row.high_seed=seed;row.phase_index=rep;row.phase_injected_deg=phase*180/pi;
row.eo_est=fit.eo_id;row.eo_success=(fit.eo_id==10);row.g_est_mm=fit.g_used;
row.A_est_mm=fit.A_id;row.rmse_V=fit.rmse;
end

function row=empty_row()
row=struct('case_id',NaN,'axis',"",'g_truth_mm',NaN,'A_true_mm',NaN,'high_seed',NaN,...
    'phase_index',NaN,'phase_injected_deg',NaN,'eo_est',NaN,'eo_success',false,...
    'g_est_mm',NaN,'A_est_mm',NaN,'rmse_V',NaN);
end
