function cfg = resolve_single_sync_experiment_paths_20250527(cfg)
%RESOLVE_SINGLE_SYNC_EXPERIMENT_PATHS_20250527 Resolve output paths.

sensor_tag = format_sensor_tag_local(cfg.analysis_sensors);
time_tag = format_time_tag_local(cfg.analysis_start_time);
case_tag = regexprep(cfg.dynamic_case_name, '[^0-9A-Za-z]+', '_');

cfg.case_output_dir = fullfile(cfg.base_cfg.output_root, cfg.dynamic_case_name);
cfg.reference_output_dir = cfg.base_cfg.reference_output_dir;
cfg.case_file = fullfile(cfg.project_dir, ...
    sprintf('ExperimentCase_20250527_B%d_%s_%s_Start%ss.mat', ...
    cfg.target_blade, sensor_tag, case_tag, time_tag));

if isfield(cfg, 'custom_calib_file') && ~isempty(cfg.custom_calib_file)
    cfg.calib_file = cfg.custom_calib_file;
else
    cfg.calib_file = fullfile(cfg.project_dir, ...
        sprintf('SGCalib_20250527_B%d_%s.mat', cfg.target_blade, sensor_tag));
end

cfg.legacy_calib_file = fullfile(cfg.project_dir, ...
    sprintf('Blade_%d_SG_Static_Config_20250527.mat', cfg.target_blade));
cfg.calib_summary_csv = fullfile(cfg.project_dir, ...
    sprintf('SGCalib_20250527_B%d_%s_Summary.csv', cfg.target_blade, sensor_tag));
cfg.result_file = fullfile(cfg.project_dir, ...
    sprintf('Result_single_sync_20250527_B%d_%s_%s_Start%ss.mat', ...
    cfg.target_blade, sensor_tag, case_tag, time_tag));
cfg.report_file = fullfile(cfg.project_dir, ...
    sprintf('Result_single_sync_20250527_B%d_%s_%s_Start%ss.txt', ...
    cfg.target_blade, sensor_tag, case_tag, time_tag));
end


function tag = format_sensor_tag_local(sensor_ids)
sensor_ids = sensor_ids(:).';
parts = compose('%d', sensor_ids);
tag = ['S', strjoin(cellstr(parts), '')];
end


function tag = format_time_tag_local(t_sec)
tag = strrep(sprintf('%.1f', t_sec), '.', 'p');
end
