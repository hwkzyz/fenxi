function Test_R5_NoDynamicStaticTilt()
%TEST_R5_NODYNAMICSTATICTILT Verify frozen-gap dynamic behavior in both cases.
repo = fileparts(mfilename('fullpath'));
common = fullfile(repo,'common');
case27 = fullfile(repo,'20250527_gap_vibration_identification_v1','latest_programs_20260907');
case22 = fullfile(repo,'20251222_gap_vibration_identification_v1','latest_programs_20260905');
addpath(common,'-begin');

% Both formal configurations must disable dynamic reconstruction from mu/tau.
addpath(case27,'-begin');
c27 = Config_20250527();
assert(isequal(double(c27.frequency.searchHz),[300 1000]), ...
    'R5:FrequencyContract','20250527 search band changed.');
assert(~logical(c27.r5.useStaticGapSlopeInDynamic), ...
    'R5:DynamicTilt','20250527 enables dynamic static-slope reconstruction.');
assert(strcmp(c27.r5.staticTiltCarrier,'absolute_surface_state'), ...
    'R5:StaticTiltCarrier','20250527 does not use the frozen absolute surface.');
m27 = R5_Build_SensorConditionedDynamicModels_20250527( ...
    fullfile(c27.paths.results,'r5_sensor_conditioned_sidecar.mat'), ...
    c27.files.lowSpeedTemplate,c27.case.targetBlade);
check_model_behavior(m27,'20250527');
rmpath(case27);
clear Config_20250527 R5_Build_SensorConditionedDynamicModels_20250527

addpath(case22,'-begin');
c22 = Config_20251222();
assert(isequal(double(c22.frequency.searchHz),[300 1000]), ...
    'R5:FrequencyContract','20251222 search band changed.');
assert(~logical(c22.r5.useStaticGapSlopeInDynamic), ...
    'R5:DynamicTilt','20251222 enables dynamic static-slope reconstruction.');
assert(strcmp(c22.r5.staticTiltCarrier,'absolute_surface_state'), ...
    'R5:StaticTiltCarrier','20251222 does not use the frozen absolute surface.');
sidecar22 = fullfile(c22.paths.results,'R5_SensorConditionedSidecar_B1_S123.mat');
if exist(sidecar22,'file') ~= 2
    error('R5:MissingSidecar','20251222 formal sidecar is missing: %s',sidecar22);
end
try
    m22 = R5_Build_SensorConditionedDynamicModels_20251222( ...
        sidecar22,c22.files.lowSpeedTemplate,c22.case.targetBlade);
catch ME
    error('R5:SidecarContractNotRebuilt', ...
        ['20251222 sidecar does not satisfy the new frozen-surface contract. ' ...
         'Rebuild it from its original R4 source before dynamic fitting.\n%s'],ME.message);
end
check_model_behavior(m22,'20251222');
fprintf('Test_R5_NoDynamicStaticTilt PASS\n');
end

function check_model_behavior(models,label)
assert(~isempty(models),'R5:NoModels','%s produced no models.',label);
for k = 1:numel(models)
    x = linspace(models(k).templateDomainMm(1), ...
        models(k).templateDomainMm(2),5).';
    dg = 0.01;
    [~,a] = models(k).evaluate(dg,x);
    [~,b] = models(k).evaluate(2*dg,x);
    assert(all(isfinite(a.gapLowMm)) && all(isfinite(a.gapHighMm)), ...
        'R5:NonFiniteGap','%s S%d returned non-finite gaps.',label,models(k).sensorId);
    assert(max(abs((a.gapHighMm-a.gapLowMm)-dg)) < 1e-12, ...
        'R5:GapIncrement','%s S%d does not use scalar dg.',label,models(k).sensorId);
    assert(max(abs((b.gapHighMm-b.gapLowMm)-2*dg)) < 1e-12, ...
        'R5:GapIncrement','%s S%d does not scale with dg.',label,models(k).sensorId);
end
end
