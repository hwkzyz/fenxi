function Diagnostic_Run_FinalContract_20250527()
% Run the matched strict-Top3/core-only Foundation -> GapAware diagnostic chain.

packageRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
foundationFile = fullfile(packageRoot, 'results', 'fixed_gap', ...
    'diag_contract_top3_coreonly', '20250526_2500-3500_t400', ...
    'FixedGap_B1_S136_diag_contract_top3_coreonly.mat');
assert(isfile(foundationFile), 'Missing diagnostic Foundation result: %s', foundationFile);

setenv('STEP07J_FOUNDATION_STEP05_RESULT_FILE', foundationFile);
setenv('STEP07J_PHASE_SAFE_EXPANSION', '0');
setenv('STEP07J_PREV_WINDOW_CANDIDATE_MODE', 'off');
setenv('STEP07J_RESULT_SUFFIX', '_diag_finalcontract_top3_coreonly');
setenv('STEP07J_MAX_WINDOWS', 'inf');
setenv('STEP07J_RUN_MODE', 'main');
run(fullfile(packageRoot, 'Main10_GapAware_FullWave_Identification_20250527.m'));
end
