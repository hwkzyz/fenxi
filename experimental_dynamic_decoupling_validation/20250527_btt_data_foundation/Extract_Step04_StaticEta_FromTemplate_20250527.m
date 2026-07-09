clc; close all;

%EXTRACT_STEP04_STATICETA_FROMTEMPLATE_20250527
% Directly compute static eta_s from Step04 template xc metadata.
%
% Formal definition:
%   eta_s = xc(reference sensor) - xc(sensor s)
%
% This is the preferred static-eta source.  Low-speed sliding-window eta is
% kept only as an independent diagnostic, not as the formal prior.

cfg = BTTProjectConfig_20250527();

S = struct();
S.blade_id = 1;
S.sensor_ids = cfg.sensor_ids;
S.reference_sensor_id = S.sensor_ids(1);
S.output_dir = fullfile(cfg.output_root, 'step05_static_eta_from_step04_xc');

if exist(S.output_dir, 'dir') ~= 7
    mkdir(S.output_dir);
end

template_file = fullfile(cfg.step04_output_dir, ...
    sprintf('Template_OPRCenterStd_LowSpeed_AllBlades_S%s_%s.mat', ...
    sprintf('%d', cfg.sensor_ids), cfg.dataset));

if ~isfile(template_file)
    error('Missing Step04 template: %s', template_file);
end

loaded = load(template_file, 'Template');
Template = loaded.Template;

[StaticEtaTable, eta_vec] = build_step04_static_eta_table_local(Template, S);

csv_file = fullfile(S.output_dir, ...
    sprintf('Step04_StaticEta_FromTemplateXc_B%d_S%s_%s.csv', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset));
mat_file = fullfile(S.output_dir, ...
    sprintf('Step04_StaticEta_FromTemplateXc_B%d_S%s_%s.mat', ...
    S.blade_id, sprintf('%d', S.sensor_ids), cfg.dataset));

writetable(StaticEtaTable, csv_file);
save(mat_file, 'StaticEtaTable', 'eta_vec', 'S', 'cfg', 'Template', ...
    'template_file', '-v7.3');

fprintf('\n=== Step04 template-xc static eta ===\n');
fprintf('Dataset: %s, blade B%d, sensors S%s, reference CH%d\n', ...
    cfg.dataset, S.blade_id, sprintf('%d', S.sensor_ids), ...
    S.reference_sensor_id);
fprintf('eta_s = xc(CH%d) - xc(CHs): %s mm\n', ...
    S.reference_sensor_id, mat2str(eta_vec, 8));
fprintf('Saved table:\n  %s\n', csv_file);
disp(StaticEtaTable);


function [T, eta_vec] = build_step04_static_eta_table_local(Template, S)
sensors = S.sensor_ids(:).';
centers = nan(numel(sensors), 1);
source = strings(numel(sensors), 1);

for i = 1:numel(sensors)
    tpl = get_template_entry_local(Template, sensors(i), S.blade_id);
    if isempty(tpl)
        error('Missing Step04 template entry for CH%d B%d.', ...
            sensors(i), S.blade_id);
    end

    [centers(i), source(i)] = get_template_xc_local(tpl);
    if ~isfinite(centers(i))
        error('Invalid Step04 xc for CH%d B%d.', sensors(i), S.blade_id);
    end
end

ref_idx = find(sensors == S.reference_sensor_id, 1, 'first');
if isempty(ref_idx)
    error('Reference sensor CH%d is not in sensor list.', S.reference_sensor_id);
end

eta_vec = centers(ref_idx).' - centers.';
eta_vec(ref_idx) = 0;

T = table(sensors(:), repmat(S.reference_sensor_id, numel(sensors), 1), ...
    repmat(S.blade_id, numel(sensors), 1), centers, source, eta_vec(:), ...
    'VariableNames', {'sensor_id', 'reference_sensor_id', 'blade_id', ...
    'template_xc_mm', 'xc_source', 'eta_prior_mm'});
end


function tpl = get_template_entry_local(Template, sid, blade_id)
tpl = [];
if ~isfield(Template, 'SensorBlade')
    return;
end

for i = 1:numel(Template.SensorBlade)
    if Template.SensorBlade(i).sensor_id == sid && ...
            Template.SensorBlade(i).blade_id == blade_id
        tpl = Template.SensorBlade(i);
        return;
    end
end
end


function [xc, source] = get_template_xc_local(tpl)
if isfield(tpl, 'xc') && isscalar(tpl.xc) && isfinite(tpl.xc)
    xc = tpl.xc;
    source = "SensorBlade.xc";
    return;
end

if isfield(tpl, 'xc_mm') && isscalar(tpl.xc_mm) && isfinite(tpl.xc_mm)
    xc = tpl.xc_mm;
    source = "SensorBlade.xc_mm";
    return;
end

if isfield(tpl, 'xc_abs') && isscalar(tpl.xc_abs) && isfinite(tpl.xc_abs)
    xc = tpl.xc_abs;
    source = "SensorBlade.xc_abs";
    return;
end

[xc, source] = estimate_template_xc_fallback_local(tpl);
end


function [xc, source] = estimate_template_xc_fallback_local(tpl)
xc = NaN;
source = "fallback_template_top85_weighted_centroid_xc";

if ~isfield(tpl, 'x_grid') || isempty(tpl.x_grid) || ...
        ~isfield(tpl, 'v_grid') || isempty(tpl.v_grid)
    return;
end

x = tpl.x_grid(:);
v = tpl.v_grid(:);
if isfield(tpl, 'v_grid_baseline_removed') && ...
        numel(tpl.v_grid_baseline_removed) == numel(x)
    vz = tpl.v_grid_baseline_removed(:);
elseif isfield(tpl, 'baseline') && isscalar(tpl.baseline) && isfinite(tpl.baseline)
    vz = v - tpl.baseline;
else
    vz = v - median(v, 'omitnan');
end

mask = isfinite(x) & isfinite(vz);
if isfield(tpl, 'domain_effective_mask') && ...
        numel(tpl.domain_effective_mask) == numel(x)
    mask = mask & tpl.domain_effective_mask(:);
end

x = x(mask);
vz = max(vz(mask), 0);

wp = vz(isfinite(vz) & vz > 0);
if numel(wp) < 5
    return;
end

high = vz >= prctile(wp, 85);
den = sum(vz(high), 'omitnan');
if den <= eps
    xc = median(x(high), 'omitnan');
else
    xc = sum(x(high) .* vz(high), 'omitnan') ./ den;
end
end
