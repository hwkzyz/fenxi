function Diagnostic_Run_FinalContract_20251222(regionId)
% Run one matched strict-Top3/core-only Foundation -> GapAware diagnostic chain.

arguments
    regionId (1,1) double {mustBeMember(regionId, [1 4 7])}
end

packageRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
switch regionId
    case 1
        bladeId = 1;
    case 4
        bladeId = 5;
    case 7
        bladeId = 1;
end
regionTag = sprintf('r%02d', regionId);
foundationFile = fullfile(packageRoot, 'results', 'fixed_gap', ...
    ['step05_noeta_vptopk_direct_template_dtop3core_', regionTag], ...
    '1000_2500_3500', sprintf( ...
    'Result_Step05_NoEtaVPTopKDirectTemplate_B%d_S123_20251222.mat', bladeId));
assert(isfile(foundationFile), 'Missing diagnostic Foundation result: %s', foundationFile);

setenv('BLADE_RESONANCE_REGION_ID', sprintf('%d', regionId));
setenv('STEP07J_FOUNDATION_STEP05_RESULT_FILE', foundationFile);
setenv('STEP07J_PHASE_SAFE_EXPANSION', '0');
setenv('STEP07J_RESULT_SUFFIX', ['_diag_finalcontract_top3_core_', regionTag]);
setenv('STEP07J_MAX_WINDOWS', 'inf');
setenv('STEP07J_RUN_MODE', 'main');
run(fullfile(packageRoot, 'Main10_GapAware_FullWave_Identification_20251222.m'));
end
