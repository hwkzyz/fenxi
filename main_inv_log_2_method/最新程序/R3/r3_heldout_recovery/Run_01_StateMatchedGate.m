function Result = Run_01_StateMatchedGate(maxCases)
%RUN_01_STATEMATCHEDGATE Noiseless calibrated-family state-matched gate.

if nargin<1||isempty(maxCases),maxCases=Inf;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(fullfile(root,'local_func'),'-begin');P=R3_Protocol();R3=r3_build_context(mainDir,P);
outDir=fullfile(root,'output','01_state_matched_gate');if ~exist(outDir,'dir'),mkdir(outDir);end
[C0,~]=r3_freeze_reference_calibration(R3,0,P.mainLowSeed);
rows=repmat(gate_row(),numel(P.gOperatingMm)*numel(P.calibratedStateMatchedGateAmplitudeMm),1);ir=0;
for g=P.gOperatingMm
    for A=P.calibratedStateMatchedGateAmplitudeMm
        if ir>=maxCases,break;end
        ir=ir+1;phi=P.referencePhaseRad;
        [~,map,truth]=r3_generate_high_case(R3,g,A,phi,0,P.mainHighSeeds(ir));
        cfg=R3.cfg;cfg.route30GapHalfWidthMm=0;cfg.route30GapCount=1;cfg.route30SingleGapHalfWidthMm=0;
        state=C0.lowState;state.gHat=g;fit=run_single_sync_voltage_vp(map,C0.templateModel,cfg,state);
        fit.r3_method="state_matched";fit.elapsed_r3_s=NaN;
        q=r3_pack_result(fit,truth,P,NaN,C0);
        rows(ir)=struct('gap_mm',g,'A_true_mm',A,'eo_est',q.eo_est,...
            'A_est_mm',q.A_est_mm,'amplitude_error_mm',q.amplitude_error_mm,...
            'phase_error_rad',q.phase_error_rad,'rmse_V',q.rmse_V,...
            'vib_success',q.vib_success,'gap_fixed_error_mm',q.gap_error_mm);
    end
    if ir>=maxCases,break;end
end
Gate=struct2table(rows(1:ir));writetable(Gate,fullfile(outDir,'state_matched_gate.csv'));
Result=struct('P',P,'Gate',Gate,'all_nonzero_success',all(Gate.vib_success(Gate.A_true_mm>0)),...
    'outputDir',outDir);save(fullfile(outDir,'state_matched_gate.mat'),'Result','-v7.3');disp(Gate);
end
function r=gate_row(),r=struct('gap_mm',NaN,'A_true_mm',NaN,'eo_est',NaN,...
    'A_est_mm',NaN,'amplitude_error_mm',NaN,'phase_error_rad',NaN,...
    'rmse_V',NaN,'vib_success',false,'gap_fixed_error_mm',NaN);end
