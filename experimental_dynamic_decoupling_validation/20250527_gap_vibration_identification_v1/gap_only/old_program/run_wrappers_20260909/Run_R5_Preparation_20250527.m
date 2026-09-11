function manifest = Run_R5_Preparation_20250527()
%RUN_R5_PREPARATION_20250527 Build the portable R5 calibration artifacts.
cfg = Setup_Paths_20250527();
check = Check_R5_Migration_20250527();
assert(check.pass, '%s', strjoin(check.messages, newline));
r5Dir = fullfile(cfg.paths.analysis, 'r5_anchor_guided_experimental_20260903');
familyFile = fullfile(cfg.paths.results, 'r5_experimental_surface_family.mat');
registrationFile = fullfile(cfg.paths.results, 'r5_b2_anchor_registration.mat');
localizationDir = fullfile(cfg.paths.results, 'r5_surface_localization');
outputFile = fullfile(cfg.paths.results, 'r5_main10_compatible_library.mat');
addpath(r5Dir);
addpath(cfg.paths.analysis);
R5_Build_R4_ExperimentalSurfaceFamily(cfg.files.responseSurface, familyFile, cfg.case.targetBlade);
% Keep the R5 B2 registration as the active experimental-family coordinate
% registration. The legacy P2 table is imported separately for comparison;
% it was fitted against the legacy production surface and cannot be silently
% substituted into the R5 family without a cross-surface equivalence check.
R5_Calibrate_B2_AnchorRegistration(familyFile, cfg.files.lowSpeedTemplate, registrationFile, 1.0, cfg.r5.anchorBlade, cfg.case.gapSensors);
R5_Localize_R4_Surface_FromLowSpeed(familyFile, cfg.files.lowSpeedTemplate, ...
    localizationDir, 1.0, registrationFile, cfg.r5.anchorBlade, ...
    cfg.case.gapSensors, cfg.r5.latentMode);
sidecarFile = fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat');
R5_Build_SensorConditionedSidecar_20250527(familyFile,fullfile(localizationDir,'loc.mat'),registrationFile,sidecarFile);
R5_Prepare_SensorConditionedDynamicContract_20250527(sidecarFile, ...
    fullfile(cfg.paths.preparedInputs,'foundation','Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat'));
exportStatus = 'completed';
exportMessage = '';
try
    R5_Export_Main10_CompatibleLibrary(familyFile, fullfile(localizationDir, 'loc.mat'), ...
        cfg.files.gapLibrary, outputFile, cfg.case.targetBlade);
catch ME
    if strcmp(ME.identifier,'R5:SensorConditionedExportRequired')
        exportStatus = 'not_applicable_sensor_conditioned';
    else
        exportStatus = 'failed';
    end
    exportMessage = ME.message;
end

manifest = struct('dataset', cfg.dataset, 'targetBlade', cfg.case.targetBlade, ...
    'responseSurface', cfg.files.responseSurface, 'template', cfg.files.lowSpeedTemplate, ...
    'baseGapLibrary', cfg.files.gapLibrary, 'familyFile', familyFile, ...
    'registrationFile', registrationFile, 'localizationDir', localizationDir, ...
    'outputFile', outputFile, 'sidecarFile', sidecarFile, 'exportStatus', exportStatus, ...
    'exportMessage', exportMessage, 'status', exportStatus);
save(fullfile(cfg.paths.results, 'r5_preparation_manifest.mat'), 'manifest', '-v7.3');
end


