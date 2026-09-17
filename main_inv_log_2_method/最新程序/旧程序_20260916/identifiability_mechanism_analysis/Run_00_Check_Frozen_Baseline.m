function Report = Run_00_Check_Frozen_Baseline()
%RUN_00_CHECK_FROZEN_BASELINE Verify read-only analysis dependencies.
analysisDir = fileparts(mfilename('fullpath'));
programDir = fileparts(analysisDir);
required = {
    'calibrate_inv_log_2_low_speed.m'
    'run_inv_log_2_structured_main.m'
    fullfile('local_func','simulate_highspeed_from_low_increment.m')
    fullfile('local_func','build_low_speed_templates_adaptive.m')
    fullfile('local_func','run_dual_sync_voltage_vp_funnel.m')
    fullfile('low_high_gap_fast_bridge','Run_DualSyncSupportAwarePostFreezeValidation.m')
    fullfile('low_high_gap_fast_bridge','Run_X01DeterministicBiasAudit.m')
    fullfile('low_high_gap_fast_bridge','Run_X02MechanismStudy.m')};
existsFlag = false(size(required));
for i = 1:numel(required)
    existsFlag(i) = isfile(fullfile(programDir,required{i}));
end
Report = table(string(required),existsFlag,'VariableNames',{'relative_path','exists'});
disp(Report);
if ~all(existsFlag)
    error('identifiability:MissingFrozenDependency', ...
        'At least one frozen dependency is missing. No analysis was run.');
end
fprintf('Frozen dependency check passed. No source file was modified.\n');
end
