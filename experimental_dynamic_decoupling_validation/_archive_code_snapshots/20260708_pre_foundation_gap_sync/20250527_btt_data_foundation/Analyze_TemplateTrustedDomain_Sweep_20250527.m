function Analyze_TemplateTrustedDomain_Sweep_20250527(mode)
%ANALYZE_TEMPLATETRUSTEDDOMAIN_SWEEP_20250527
% Build and summarize trusted-domain variants for Step05 sensitivity tests.
%
% The first diagnostic variant keeps the foundation template curve unchanged
% but replaces its trusted x-domain with the narrower calibration domain.
% This isolates trusted-domain width from template-curve differences.

if nargin < 1 || isempty(mode)
    mode = 'prepare';
end
mode = lower(strtrim(char(mode)));

clc;
P = settings_local();
ensure_dir_local(P.out_dir);

switch mode
    case 'prepare'
        prepare_domain_variants_local(P);
    case 'summarize'
        summarize_domain_sweep_local(P);
    otherwise
        error('Unknown mode: %s. Use prepare or summarize.', mode);
end
end


function P = settings_local()
script_dir = fileparts(mfilename('fullpath'));
cfg = BTTProjectConfig_20250527();

P = struct();
P.cfg = cfg;
P.script_dir = script_dir;
P.dataset = cfg.dataset;
P.case_name = '20250526_2500-3500_t400';
P.blade_id = 1;
P.sensors = [1 3 6];
P.sensor_tag = ['S', sprintf('%d', P.sensors)];
P.out_dir = fullfile(cfg.output_root, 'template_trusted_domain_sweep', P.case_name);

P.foundation_template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_%s_%s.mat', ...
    P.sensor_tag, P.dataset));
P.calibration_adapted_template_file = fullfile(cfg.output_root, ...
    'template_calibration_impact', P.case_name, sprintf( ...
    'Template_CalibrationAdapted_ForFoundationStep05_B%d_%s_%s.mat', ...
    P.blade_id, P.sensor_tag, P.dataset));

P.variant_specs = struct( ...
    'tag', {'domF95', 'domF90', 'domF88', 'domF86'}, ...
    'mode', {'ratio', 'ratio', 'ratio', 'calibration_domain'}, ...
    'ratio', {0.95, 0.90, 0.88, NaN});
P.variant_summary_csv = fullfile(P.out_dir, sprintf( ...
    'TrustedDomainSweep_TemplateVariants_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
P.step05_summary_csv = fullfile(P.out_dir, sprintf( ...
    'TrustedDomainSweep_Step05Summary_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
P.window_delta_csv = fullfile(P.out_dir, sprintf( ...
    'TrustedDomainSweep_WindowDelta_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));

P.baseline_step05_root = fullfile(cfg.output_root, ...
    'step05_single_sync_direct_template_impF');
end


function prepare_domain_variants_local(P)
if ~isfile(P.foundation_template_file)
    error('Foundation template not found: %s', P.foundation_template_file);
end
if ~isfile(P.calibration_adapted_template_file)
    error(['Calibration-adapted template not found: %s\n' ...
        'Run Analyze_TemplateCalibration_Impact_20250527(''prepare'') first.'], ...
        P.calibration_adapted_template_file);
end

rows = {};
F = load(P.foundation_template_file, 'Template');
C = load(P.calibration_adapted_template_file, 'Template');

for iv = 1:numel(P.variant_specs)
    spec = P.variant_specs(iv);
    Variant = F.Template;
    Variant.CreatedBy = mfilename;
    Variant.CreatedOn = datestr(now, 31);
    Variant.TrustedDomainSweep_Tag = spec.tag;
    Variant.TrustedDomainSweep_Mode = spec.mode;

    for i = 1:numel(P.sensors)
        sid = P.sensors(i);
        idxF = find_entry_index_local(Variant, sid, P.blade_id);
        idxC = find_entry_index_local(C.Template, sid, P.blade_id);
        if isempty(idxF) || isempty(idxC)
            error('Missing template entry CH%d B%d.', sid, P.blade_id);
        end

        eF = Variant.SensorBlade(idxF);
        eC = C.Template.SensorBlade(idxC);
        old_domain = eF.x_domain;

        switch lower(spec.mode)
            case 'ratio'
                center = mean(old_domain);
                half_width = 0.5 * diff(old_domain) * spec.ratio;
                new_domain = [center - half_width, center + half_width];
                quality = sprintf('foundation_curve_domain_ratio_%.3f', spec.ratio);
            case 'calibration_domain'
                new_domain = eC.x_domain;
                quality = 'foundation_curve_calibration_domain';
            otherwise
                error('Unknown domain variant mode: %s', spec.mode);
        end

        eF.x_domain = new_domain;
        eF.domain_mask = eF.x_grid(:) >= new_domain(1) & ...
            eF.x_grid(:) <= new_domain(2) & ...
            isfinite(eF.x_grid(:)) & isfinite(eF.v_grid(:));
        eF.domain_effective_mask = eF.domain_mask;
        eF.quality_status = quality;
        Variant.SensorBlade(idxF) = eF;

        rows(end+1, :) = {string(spec.tag), string(spec.mode), spec.ratio, ...
            sid, P.blade_id, old_domain(1), old_domain(2), ...
            diff(old_domain), new_domain(1), new_domain(2), ...
            diff(new_domain), diff(new_domain) / diff(old_domain)}; %#ok<AGROW>
    end

    Template = Variant; %#ok<NASGU>
    metadata = struct();
    metadata.created_by = mfilename;
    metadata.created_on = datestr(now, 31);
    metadata.variant_tag = spec.tag;
    metadata.variant_mode = spec.mode;
    metadata.variant_ratio = spec.ratio;
    variant_file = variant_template_file_local(P, spec.tag);
    save(variant_file, 'Template', 'metadata', '-v7.3');
    fprintf('Trusted-domain variant saved [%s]:\n  %s\n', ...
        spec.tag, variant_file);
end

T = cell2table(rows, 'VariableNames', { ...
    'variant_tag', 'variant_mode', 'variant_ratio', ...
    'sensor_id', 'blade_id', ...
    'foundation_left_mm', 'foundation_right_mm', 'foundation_width_mm', ...
    'variant_left_mm', 'variant_right_mm', 'variant_width_mm', ...
    'variant_width_ratio'});
writetable(T, P.variant_summary_csv);

fprintf('Variant summary saved:\n  %s\n', P.variant_summary_csv);
fprintf('\nRun Step05 tests with STEP05_OUTPUT_TAG equal to domF95/domF90/domF86.\n');
end


function summarize_domain_sweep_local(P)
baseline_file = step05_trend_file_local(P, P.baseline_step05_root);
summary_rows = {};
all_delta = table();

if ~isfile(baseline_file)
    error('Baseline Step05 file not found: %s', baseline_file);
end

B = readtable(baseline_file);

for iv = 1:numel(P.variant_specs)
    tag = P.variant_specs(iv).tag;
    root = fullfile(P.cfg.output_root, ['step05_single_sync_direct_template_', tag]);
    variant_file = step05_trend_file_local(P, root);
    if ~isfile(variant_file)
        fprintf('Variant Step05 file missing, skip %s:\n  %s\n', tag, variant_file);
        continue;
    end

    V = readtable(variant_file);
    D = build_window_delta_local(B, V, tag);
    all_delta = [all_delta; D]; %#ok<AGROW>

    summary_rows(end+1, :) = {string(tag), height(D), nnz(D.EO_changed), ...
        mean(D.A_baseline_mm, 'omitnan'), ...
        mean(D.A_variant_mm, 'omitnan'), ...
        mean(D.A_delta_mm, 'omitnan'), ...
        mean(D.fn_baseline_hz, 'omitnan'), ...
        mean(D.fn_variant_hz, 'omitnan'), ...
        mean(D.rmse_baseline_mv, 'omitnan'), ...
        mean(D.rmse_variant_mv, 'omitnan'), ...
        mean(D.rmse_delta_mv, 'omitnan'), ...
        mean(D.point_count_ratio, 'omitnan'), ...
        mean(D.core_point_count_ratio, 'omitnan')}; %#ok<AGROW>
end

if isempty(summary_rows)
    error('No variant Step05 outputs found.');
end

writetable(all_delta, P.window_delta_csv);

S = cell2table(summary_rows, 'VariableNames', { ...
    'variant_tag', 'n_windows', 'n_eo_changed', ...
    'mean_A_baseline_mm', 'mean_A_variant_mm', 'mean_A_delta_mm', ...
    'mean_fn_baseline_hz', 'mean_fn_variant_hz', ...
    'mean_rmse_baseline_mv', 'mean_rmse_variant_mv', ...
    'mean_rmse_delta_mv', 'mean_point_count_ratio', ...
    'mean_core_point_count_ratio'});
writetable(S, P.step05_summary_csv);

fprintf('Trusted-domain Step05 delta saved:\n  %s\n', P.window_delta_csv);
fprintf('Trusted-domain Step05 summary saved:\n  %s\n', P.step05_summary_csv);
disp(S);
end


function D = build_window_delta_local(B, V, tag)
common = intersect(B.window_id, V.window_id, 'stable');
n = numel(common);

D = table();
D.variant_tag = repmat(string(tag), n, 1);
D.window_id = common(:);
D.EO_baseline = nan(n,1);
D.EO_variant = nan(n,1);
D.A_baseline_mm = nan(n,1);
D.A_variant_mm = nan(n,1);
D.fn_baseline_hz = nan(n,1);
D.fn_variant_hz = nan(n,1);
D.rmse_baseline_mv = nan(n,1);
D.rmse_variant_mv = nan(n,1);
D.point_count_baseline = nan(n,1);
D.point_count_variant = nan(n,1);
D.core_point_count_baseline = nan(n,1);
D.core_point_count_variant = nan(n,1);

for i = 1:n
    wb = B(B.window_id == common(i), :);
    wv = V(V.window_id == common(i), :);
    D.EO_baseline(i) = value_local(wb, 'EO_id');
    D.EO_variant(i) = value_local(wv, 'EO_id');
    D.A_baseline_mm(i) = abs(value_local(wb, 'A_id'));
    D.A_variant_mm(i) = abs(value_local(wv, 'A_id'));
    D.fn_baseline_hz(i) = value_local(wb, 'fn_id');
    D.fn_variant_hz(i) = value_local(wv, 'fn_id');
    D.rmse_baseline_mv(i) = 1000 * value_local(wb, 'weighted_voltage_rmse');
    D.rmse_variant_mv(i) = 1000 * value_local(wv, 'weighted_voltage_rmse');
    D.point_count_baseline(i) = value_local(wb, 'point_count');
    D.point_count_variant(i) = value_local(wv, 'point_count');
    D.core_point_count_baseline(i) = value_local(wb, 'core_point_count');
    D.core_point_count_variant(i) = value_local(wv, 'core_point_count');
end

D.EO_changed = D.EO_variant ~= D.EO_baseline;
D.A_delta_mm = D.A_variant_mm - D.A_baseline_mm;
D.A_ratio_variant_over_baseline = D.A_variant_mm ./ D.A_baseline_mm;
D.fn_delta_hz = D.fn_variant_hz - D.fn_baseline_hz;
D.rmse_delta_mv = D.rmse_variant_mv - D.rmse_baseline_mv;
D.point_count_ratio = D.point_count_variant ./ D.point_count_baseline;
D.core_point_count_ratio = D.core_point_count_variant ./ D.core_point_count_baseline;
end


function file = variant_template_file_local(P, tag)
file = fullfile(P.out_dir, sprintf( ...
    'Template_FoundationCurve_%s_B%d_%s_%s.mat', ...
    tag, P.blade_id, P.sensor_tag, P.dataset));
end


function file = step05_trend_file_local(P, root)
file = fullfile(root, P.case_name, sprintf( ...
    'Trend_Step05_SingleSyncDirectTemplate_B%d_%s_%s.csv', ...
    P.blade_id, P.sensor_tag, P.dataset));
end


function idx = find_entry_index_local(Template, sid, blade_id)
idx = [];
for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && ...
            Template.SensorBlade(i).blade_id == blade_id
        idx = i;
        return;
    end
end
end


function v = value_local(T, name)
if ismember(name, T.Properties.VariableNames) && ~isempty(T.(name))
    v = T.(name)(1);
else
    v = NaN;
end
end


function ensure_dir_local(d)
if exist(d, 'dir') ~= 7
    mkdir(d);
end
end
