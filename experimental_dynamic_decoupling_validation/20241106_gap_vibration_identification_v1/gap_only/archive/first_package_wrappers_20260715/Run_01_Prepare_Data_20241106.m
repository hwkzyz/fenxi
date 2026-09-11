function preparedFile = Run_01_Prepare_Data_20241106()
%RUN_01_PREPARE_DATA_20241106 Build the common folder-local input layer.
cfg = Config_20241106();
add_package_paths_local(cfg);
Check_Program_20241106('code', false);
Check_Program_20241106('inputs', false);
ensure_result_dirs_local(cfg);

if cfg.preparation.rebuildFromRaw
    run_core_script_local(fullfile(cfg.paths.preparation, ...
        'Step01_Build_LowSpeed_Reference_20241106.m'));
    run_core_script_local(fullfile(cfg.paths.preparation, ...
        'Step02_Extract_Dynamic_BTT_20241106.m'));
    run_core_script_local(fullfile(cfg.paths.foundation, ...
        'Step03_Build_Dynamic_Observation_Bundle_20241106.m'));
    run_core_script_local(fullfile(cfg.paths.preparation, ...
        'Step04_Build_LowSpeed_OPRCenterStd_Template_20241106.m'));
    run_core_script_local(fullfile(cfg.paths.preparation, ...
        'Step02_Calc_BTT_Displacement_20241106.m'));
else
    copy_tree_contents_local(cfg.paths.bundledPreparedFoundation, ...
        cfg.paths.preparedFoundation);
    copy_tree_contents_local(fullfile(cfg.paths.calibrationFoundation, ...
        'step04_low_speed_template'), fullfile(cfg.paths.preparedFoundation, ...
        'step04_low_speed_template'));
    copyfile(fullfile(cfg.paths.bundledPreparedGap, ...
        'Step02_BTT_Displacement_20241106.mat'), cfg.paths.gapRuntime, 'f');
end

Check_Program_20241106('prepared', false);
caseTag = build_case_tag_local(cfg);
preparedFile = fullfile(cfg.paths.prepared, ...
    ['PreparedCase_', caseTag, '.mat']);
PreparedCase = struct(); %#ok<NASGU>
PreparedCase.packageVersion = cfg.packageVersion;
PreparedCase.createdAt = datetime('now');
PreparedCase.cfg = cfg;
PreparedCase.caseTag = caseTag;
PreparedCase.expectedWindowCount = cfg.window.count;
PreparedCase.foundationStep02File = fullfile(cfg.paths.preparedFoundation, ...
    'step02_dynamic_btt', cfg.case.dynamicCase, ...
    'Step02_Dynamic_BTT_Extraction_20241106.mat');
PreparedCase.gapStep02File = fullfile(cfg.paths.gapRuntime, ...
    'Step02_BTT_Displacement_20241106.mat');
PreparedCase.foundationTemplateFile = fullfile(cfg.paths.preparedFoundation, ...
    'step04_low_speed_template', ...
    'Template_OPRCenterStd_LowSpeed_AllBlades_S2357_20241106.mat');
PreparedCase.gapTemplateBankFile = fullfile(cfg.paths.calibrationGap, ...
    'LowSpeedTemplateBank_20241106_B1toB6_S2357.mat');
PreparedCase.gapCalibrationBankFile = fullfile(cfg.paths.calibrationGap, ...
    'GapCalibrationBank_20241106_B1toB6_S57.mat');
save(preparedFile, 'PreparedCase', '-v7.3');
fprintf('Prepared input complete:\n  %s\n', preparedFile);
end

function add_package_paths_local(cfg)
addpath(cfg.paths.root, cfg.paths.foundation, cfg.paths.preparation, ...
    cfg.paths.gapAware, cfg.paths.utilities);
end

function ensure_result_dirs_local(cfg)
dirs = {cfg.paths.prepared, cfg.paths.preparedFoundation, cfg.paths.gapRuntime, ...
    cfg.paths.foundationResults, cfg.paths.gapResults, cfg.paths.comparison, cfg.paths.figures};
for i = 1:numel(dirs)
    if ~isfolder(dirs{i})
        mkdir(dirs{i});
    end
end
end

function copy_tree_contents_local(sourceDir, targetDir)
assert(isfolder(sourceDir), 'Missing bundled prepared directory: %s', sourceDir);
if ~isfolder(targetDir)
    mkdir(targetDir);
end
items = dir(sourceDir);
items = items(~ismember({items.name}, {'.','..'}));
for i = 1:numel(items)
    sourceItem = fullfile(items(i).folder, items(i).name);
    targetItem = fullfile(targetDir, items(i).name);
    copyfile(sourceItem, targetItem, 'f');
end
end

function run_core_script_local(scriptFile)
assert(isfile(scriptFile), 'Missing core script: %s', scriptFile);
run(scriptFile);
end

function tag = build_case_tag_local(cfg)
sensorTag = ['S', sprintf('%d', cfg.case.analysisSensors)];
timeTag = strrep(sprintf('T%07.3f', cfg.case.analysisStartTimeSec), '.', 'p');
tag = sprintf('B%d_%s_%s_W%dS%d_F%dto%d', cfg.case.targetBlade, sensorTag, ...
    timeTag, cfg.window.windowBladePasses, cfg.window.slidingStepBladePasses, ...
    round(cfg.frequency.searchHz(1)), round(cfg.frequency.searchHz(2)));
end
