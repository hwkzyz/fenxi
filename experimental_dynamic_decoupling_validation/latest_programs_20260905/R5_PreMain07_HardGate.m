function Gate = R5_PreMain07_HardGate(resultFile)
%R5_PREMAIN07_HARDGATE Hard provenance/data gate before strain validation.
% Scientific quality concerns remain flags; only integrity failures block.
if nargin < 1 || isempty(resultFile)
    error('R5:GateInput','Pass the newly generated Main06 MAT result explicitly.');
end
S = load(resultFile,'Result');
assert(isfield(S,'Result'),'R5:GateInput','MAT file lacks Result.');
R = S.Result; reasons = strings(0,1); flags = strings(0,1);
if ~isfield(R,'MethodContract'), reasons(end+1)="missing_method_contract"; end %#ok<AGROW>
formalMethod = 'R5_PC1_AnchorSurface_SensorWiseDg_FullWave';
if ~isfield(R,'method') || ~strcmp(string(R.method),formalMethod)
    reasons(end+1)="wrong_formal_method"; %#ok<AGROW>
end
if ~isfield(R,'cfg') || ~isfield(R,'flowConfig') || ...
        ~isfield(R.cfg,'targetBlade') || ~isfield(R.flowConfig,'identification') || ...
        ~isfield(R.flowConfig.identification,'targetBlade') || ...
        double(R.cfg.targetBlade) ~= double(R.flowConfig.identification.targetBlade)
    reasons(end+1)="target_blade_region_mismatch"; %#ok<AGROW>
end
if ~isfield(R,'foundationStep05File') || exist(R.foundationStep05File,'file') ~= 2
    reasons(end+1)="missing_foundation_provenance"; %#ok<AGROW>
else
    F = load(R.foundationStep05File,'Result');
    if ~isfield(F,'Result') || ~isfield(F.Result,'TargetBlade') || ...
            ~isfield(R,'cfg') || double(F.Result.TargetBlade) ~= double(R.cfg.targetBlade)
        reasons(end+1)="foundation_blade_mismatch"; %#ok<AGROW>
    elseif ~isfield(F.Result,'AnalysisSettings') || ~isfield(R,'flowConfig') || ...
            abs(double(F.Result.AnalysisSettings.analysis_start_time_s) - ...
                double(R.flowConfig.identification.analysisStartTimeSec)) > 1e-9 || ...
            double(F.Result.AnalysisSettings.target_laps) ~= ...
                double(R.flowConfig.identification.targetBladePasses)
        reasons(end+1)="foundation_region_window_mismatch"; %#ok<AGROW>
    end
end
if isfield(R,'MethodContract')
    if ~isfield(R.MethodContract,'formal_method_name') || ...
            ~strcmp(string(R.MethodContract.formal_method_name),formalMethod)
        reasons(end+1)="wrong_public_method_name"; %#ok<AGROW>
    end
    if ~isfield(R.MethodContract,'front_end') || ~contains(string(R.MethodContract.front_end),'PC1')
        reasons(end+1)="wrong_front_end"; %#ok<AGROW>
    end
    if ~isfield(R.MethodContract,'sensor_wise_dynamic_gap') || ~R.MethodContract.sensor_wise_dynamic_gap
        reasons(end+1)="not_sensor_wise_dg"; %#ok<AGROW>
    end
    if isfield(R.MethodContract,'pc1_dynamic_state') && R.MethodContract.pc1_dynamic_state
        reasons(end+1)="dynamic_pc1_state_forbidden"; %#ok<AGROW>
    end
    if isfield(R.MethodContract,'delta_mu_enabled') && R.MethodContract.delta_mu_enabled
        reasons(end+1)="delta_mu_enabled"; %#ok<AGROW>
    end
    if isfield(R.MethodContract,'delta_tau_enabled') && R.MethodContract.delta_tau_enabled
        reasons(end+1)="delta_tau_enabled"; %#ok<AGROW>
    end
    if ~isfield(R.MethodContract,'data_residual') || ...
            ~strcmp(string(R.MethodContract.data_residual),'plain_unweighted_voltage_SSE') || ...
            ~isfield(R.MethodContract,'optimization_objective')
        reasons(end+1)="incomplete_objective_contract"; %#ok<AGROW>
    end
end
if ~isfield(R,'cfg') || ~isfield(R.cfg,'runMode') || ~strcmpi(string(R.cfg.runMode),'main')
    reasons(end+1)="not_formal_main_run"; %#ok<AGROW>
end
if ~isfield(R,'MainModel') || ~strcmpi(string(R.MainModel),'gap_only')
    reasons(end+1)="wrong_internal_main_branch"; %#ok<AGROW>
end
if ~isfield(R,'FoundationStep05SourceMode') || ...
        ~strcmpi(string(R.FoundationStep05SourceMode),'foundation_step05_bundle')
    reasons(end+1)="wrong_foundation_source_mode"; %#ok<AGROW>
end
if ~isfield(R,'Main06WindowLimit') || ~isinf(double(R.Main06WindowLimit))
    reasons(end+1)="limited_main06_run"; %#ok<AGROW>
end
if ~isfield(R,'FoundationWindowCoverage') || numel(R.FoundationWindowCoverage) ~= 2 || ...
        any(~isfinite(double(R.FoundationWindowCoverage))) || ...
        double(R.FoundationWindowCoverage(1)) ~= double(R.FoundationWindowCoverage(2))
    reasons(end+1)="incomplete_foundation_window_coverage"; %#ok<AGROW>
end
expectedWindows = NaN;
if isfield(R,'flowConfig') && isfield(R.flowConfig,'identification')
    I = R.flowConfig.identification;
    needed = {'targetBladePasses','windowBladePasses','slidingStepBladePasses'};
    if all(isfield(I,needed)) && I.slidingStepBladePasses > 0
        expectedWindows = floor((double(I.targetBladePasses)-double(I.windowBladePasses)) / ...
            double(I.slidingStepBladePasses)) + 1;
    end
end
if ~isfinite(expectedWindows) || ~isfield(R,'WindowResult') || ...
        numel(R.WindowResult) ~= expectedWindows || ...
        (isfield(R,'Trend') && height(R.Trend) ~= expectedWindows)
    reasons(end+1)="incomplete_main06_window_coverage"; %#ok<AGROW>
end
if ~isfield(R,'cfg') || ~isfield(R.cfg,'gapActiveSensorMask') || ...
        ~all(logical(R.cfg.gapActiveSensorMask(:))) || ...
        ~isfield(R.cfg,'analysisSensors') || numel(R.cfg.gapActiveSensorMask) ~= numel(R.cfg.analysisSensors)
    reasons(end+1)="inactive_formal_gap_sensor"; %#ok<AGROW>
end
if ~isfield(R,'ZeroGapEquivalenceGlobalMaxAbsMv') || ...
        ~isfinite(double(R.ZeroGapEquivalenceGlobalMaxAbsMv)) || ...
        abs(double(R.ZeroGapEquivalenceGlobalMaxAbsMv)) > 1e-9
    reasons(end+1)="zero_gap_equivalence_failed"; %#ok<AGROW>
end
if ~isfield(R,'CoordinateCheck') || ~isfield(R.CoordinateCheck,'is_consistent') || ...
        ~logical(R.CoordinateCheck.is_consistent)
    reasons(end+1)="coordinate_check_failed"; %#ok<AGROW>
end
if ~isfield(R,'Trend') || isempty(R.Trend), reasons(end+1)="missing_trend"; end %#ok<AGROW>
if isfield(R,'Trend')
    T = R.Trend;
    if ~ismember('gap_rmse_mV',T.Properties.VariableNames) || any(~isfinite(T.gap_rmse_mV))
        reasons(end+1)="nonfinite_rmse"; %#ok<AGROW>
    end
    if ismember('gap_EO',T.Properties.VariableNames)
        eo = double(T.gap_EO); eo = eo(isfinite(eo));
        if numel(unique(eo)) > 1, flags(end+1)="EO_branch_switch"; end %#ok<AGROW>
    end
    if ismember('gap_frequency_hz',T.Properties.VariableNames)
        f = double(T.gap_frequency_hz); f = f(isfinite(f));
        if numel(f)>1 && std(f) > 5, flags(end+1)="frequency_instability"; end %#ok<AGROW>
    end
end
if ~isfield(R,'WindowResult') || isempty(R.WindowResult)
    reasons(end+1)="missing_window_results"; %#ok<AGROW>
else
    for iw=1:numel(R.WindowResult)
        W=R.WindowResult(iw);
        if ~isfield(W,'modelFits') || ~isfield(W.modelFits,'gap_only') || ...
                ~isfield(W.modelFits.gap_only,'deltaGapMm') || ...
                any(~isfinite(W.modelFits.gap_only.deltaGapMm(:)))
            reasons(end+1)="invalid_sensorwise_dg"; break; %#ok<AGROW>
        end
        fit = W.modelFits.gap_only;
        requiredFit = {'EO','freqHz','amplitudeMm','phaseRad','dxMm','plainRmseMv','score'};
        invalidFit = ~all(isfield(fit,requiredFit));
        if ~invalidFit
            values = [fit.EO,fit.freqHz,fit.amplitudeMm,fit.phaseRad,fit.dxMm,fit.plainRmseMv,fit.score];
            invalidFit = any(~isfinite(double(values)));
        end
        if invalidFit
            reasons(end+1)="invalid_solver_or_domain_state"; break; %#ok<AGROW>
        end
    end
end
status = ternary_local(~isempty(reasons),'BLOCKED',ternary_local(~isempty(flags),'PASS_WITH_FLAGS','PASS'));
Gate = struct('status', status, 'hard_reasons', reasons, ...
    'scientific_flags', flags, 'resultFile',resultFile);
disp(Gate);
if strcmp(Gate.status,'BLOCKED'), error('R5:PreMain07Blocked','Main07 blocked: %s',strjoin(reasons,', ')); end
end
function y=ternary_local(c,a,b), if c,y=a;else,y=b;end,end
