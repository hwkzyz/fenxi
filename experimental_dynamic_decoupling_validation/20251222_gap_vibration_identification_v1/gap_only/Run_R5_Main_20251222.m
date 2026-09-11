function Run = Run_R5_Main_20251222()
%RUN_R5_MAIN_20251222 Unified R5 production entry point.
%
% The single production route is:
%   Foundation raw-window bundle -> R5 PC1/anchor sidecar ->
%   R5 PC1-Anchor SensorWiseDynamicGap Full-Wave method
%
% Main07 is deliberately not run by default.  Set R5_RUN_MAIN07=1 only after
% the returned hard gate and scientific audit have been reviewed.

rootDir = fileparts(mfilename('fullpath'));
Setup_Paths_20251222();
addpath(rootDir);
assert_production_environment_local();
[selection,~] = ResonanceRegionCatalog_20251222();
regionId = selection.regionId;
rawRegion = strtrim(getenv('BLADE_RESONANCE_REGION_ID'));
if ~isempty(rawRegion)
    regionId = str2double(rawRegion);
end
if ~isfinite(regionId) || regionId ~= round(regionId)
    error('R5:InvalidRegion','BLADE_RESONANCE_REGION_ID must be an integer resonance-region id.');
end
setenv('BLADE_RESONANCE_REGION_ID',num2str(regionId));
runMain07 = any(strcmpi(strtrim(getenv('R5_RUN_MAIN07')), {'1','true','yes','on'}));

% These are the only production stages.  Each stage reads/writes its own
% frozen contract, but the user invokes one entry point for the chain.
run(fullfile(rootDir,'Main05_Foundation_NoGap_VPTopK_20251222.m'));
run(fullfile(rootDir,'Main05_R5_PC1_FrontEnd_20251222.m'));
run(fullfile(rootDir,'Main06_GapAware_FullWave_Identification_20251222.m'));

Pflow = ProjectionFlow_Config_20251222();
targetBlade = Pflow.identification.targetBlade;
sensorTag = ['S',sprintf('%d',Pflow.identification.analysisSensors)];
resultFile = fullfile(rootDir,'results','gap_aware',sprintf( ...
    'Main_GapAware_VPTopK_FullWave_20251222_B%d_%s_%s_gapaware.mat', ...
    targetBlade,sensorTag,Pflow.identification.resonanceRegionShortTag));
if exist(resultFile,'file') ~= 2
    error('R5:MissingMain06Result','Main06 did not produce the expected result: %s',resultFile);
end

Gate = R5_PreMain07_HardGate(resultFile);
Audit = R5_Audit_SensorWiseDg(resultFile);
Comparison = R5_Compare_Proposed_vs_Foundation(resultFile);
% The stage scripts run in this function workspace and may clear ordinary
% variables. Re-resolve the frozen selection before constructing the return
% object so a successful production run cannot fail during final packaging.
[selection,~] = ResonanceRegionCatalog_20251222(regionId);
Run = struct('region',selection,'resultFile',resultFile,'gate',Gate,'audit',Audit, ...
    'comparison',Comparison, ...
    'main07_run',false);

if runMain07
    if strcmp(Gate.status,'BLOCKED')
        error('R5:Main07Blocked','Main07 is blocked by the hard gate.');
    end
    previousValidatedFile = getenv('R5_VALIDATED_RESULT_FILE');
    setenv('R5_VALIDATED_RESULT_FILE',resultFile);
    try
        run(fullfile(rootDir,'Main07_Validate_Identification_With_Strain_20251222.m'));
    catch ME
        setenv('R5_VALIDATED_RESULT_FILE',previousValidatedFile);
        rethrow(ME);
    end
    setenv('R5_VALIDATED_RESULT_FILE',previousValidatedFile);
    Run.main07_run = true;
end
end

function assert_production_environment_local()
% Reject stale smoke/diagnostic controls before any production stage runs.
% Only region selection, Main07 enablement and plot-display controls are
% permitted through the environment on the unified production route.
[status, raw] = system('set');
if status ~= 0
    error('R5:EnvironmentAuditFailed','Unable to enumerate the production environment.');
end
lines = splitlines(string(raw));
names = extractBefore(lines, '=');
names = names(strlength(names) > 0);
isStepOverride = startsWith(names,'STEP05_') | startsWith(names,'STEP07J_');
allowed = ismember(names, ["STEP05_SHOW_PLOTS","STEP05_SAVE_FIGURES"]);
forbidden = unique(names(isStepOverride & ~allowed), 'stable');
if ~isempty(strtrim(getenv('BLADE_ANALYSIS_START_TIME_SEC')))
    forbidden(end+1) = "BLADE_ANALYSIS_START_TIME_SEC"; %#ok<AGROW>
end
if ~isempty(forbidden)
    error('R5:ProductionEnvironmentOverride', ...
        ['Formal Run_R5_Main refuses numerical/path overrides left in the environment: %s. ' ...
         'Clear them and select the case only with BLADE_RESONANCE_REGION_ID.'], ...
        strjoin(forbidden, ', '));
end
end
