function manifest = Run_R5_Preparation_20241106()
%RUN_R5_PREPARATION_20241106 Build R5 calibration artifacts up to export.
cfg = Setup_Paths_20241106();
addpath(cfg.paths.analysis);
adapter = cfg.files.responseSurface;
if ~isfile(adapter)
    Build_ResponseSurface_Adapter_20241106();
end
check = Check_R5_Migration_20241106();
assert(check.pass, '%s', strjoin(check.messages, newline));
r5Dir = fullfile(cfg.paths.analysis, 'r5_anchor_guided_experimental_20260903');
addpath(r5Dir);
familyFile = fullfile(cfg.paths.results, 'r5_experimental_surface_family.mat');
registrationFile = fullfile(cfg.paths.results, 'r5_b2_anchor_registration.mat');
localizationDir = fullfile(cfg.paths.results, 'r5_surface_localization');
templateAdapter = fullfile(cfg.paths.results,'r5_low_speed_template_bank_adapter.mat');
R5_Build_LowSpeedTemplateBankAdapter_20241106([],templateAdapter);
T0 = R5_StageMatForMatlab(templateAdapter,'r5_template_adapter');
anchorBank = R5_StageMatForMatlab(cfg.files.gapLibrary,'r5_gap_bank');
ST = load(T0,'Template'); SB = load(anchorBank,'GapCalibrationBank');
anchorSensors = double(cfg.case.gapSensors(:).');
anchorGaps = nan(size(anchorSensors));
for ii=1:numel(anchorSensors)
    q = SB.GapCalibrationBank.entry([SB.GapCalibrationBank.entry.bladeId]==cfg.r5.anchorBlade & ...
        [SB.GapCalibrationBank.entry.sensorId]==anchorSensors(ii));
    assert(numel(q)==1,'R5:AnchorGapMissing','B%d/S%d anchor gap missing.',cfg.r5.anchorBlade,anchorSensors(ii));
    anchorGaps(ii)=double(q.sensor.g0Mm);
end
R5_Build_R4_ExperimentalSurfaceFamily(cfg.files.responseSurface, familyFile, cfg.case.targetBlade);
R5_Calibrate_B2_AnchorRegistration(familyFile, templateAdapter, registrationFile, anchorGaps, cfg.r5.anchorBlade, cfg.case.gapSensors);
R5_Localize_R4_Surface_FromLowSpeed(familyFile, templateAdapter, ...
    localizationDir, anchorGaps, registrationFile, cfg.r5.anchorBlade, ...
    cfg.case.gapSensors, cfg.r5.latentMode);
sidecarFile = fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat');
R5_Build_SensorConditionedSidecar_20241106(familyFile,fullfile(localizationDir,'loc.mat'),registrationFile,sidecarFile);
contract = R5_Prepare_SensorConditionedDynamicContract_20241106(sidecarFile,'',templateAdapter);
manifest = struct('dataset', cfg.dataset, 'targetBlade', cfg.case.targetBlade, ...
    'responseSurface', cfg.files.responseSurface, 'template', cfg.files.lowSpeedTemplate, ...
    'templateAdapter',templateAdapter,'anchorBlade',cfg.r5.anchorBlade,'anchorSensors',anchorSensors,'anchorGapsMm',anchorGaps, ...
    'familyFile', familyFile, 'registrationFile', registrationFile, ...
    'localizationDir', localizationDir, 'sidecarFile',sidecarFile, ...
    'dynamicContractFile',fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_contract.mat'), ...
    'status', 'sensor_conditioned_prepared');
save(fullfile(cfg.paths.results, 'r5_preparation_manifest.mat'), 'manifest', '-v7.3');
end


