%% Step01B_Build_GradientXRange030_BaseFrame_Template_20250527
% Represent the GradientXRange030 waveform in the base-template x frame.
%
% The GradientXRange030 template keeps the high-sensitivity waveform domain
% and weights.  This diagnostic/theory-clean variant shifts its x_grid and
% x_domain into the base DynamicMap coordinate frame, so Step03 interpolation
% uses the same x-zero convention as DynamicMap_B1_S136_SlidingWindows.

clear; clc; close all;

route_dir = fileparts(mfilename('fullpath'));
sensor_tag = 'S136';
date_tag = '20250527';
target_blade = 1;
analysis_sensors = [1 3 6];

template_dir = fullfile(route_dir, 'output', 'templates');
base_template_file = fullfile(template_dir, ...
    sprintf('Template_LowSpeedRotating_B%d_%s_%s.mat', target_blade, sensor_tag, date_tag));
gradient_template_file = fullfile(template_dir, ...
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_%s.mat', target_blade, sensor_tag, date_tag));
base_frame_template_file = fullfile(template_dir, ...
    sprintf('Template_LowSpeedRotating_B%d_%s_GradientXRange030_BaseFrame_%s.mat', target_blade, sensor_tag, date_tag));

if exist(base_template_file, 'file') ~= 2
    error('Base template not found: %s', base_template_file);
end
if exist(gradient_template_file, 'file') ~= 2
    error('GradientXRange030 template not found: %s', gradient_template_file);
end

loaded_base = load(base_template_file, 'Template');
loaded_gradient = load(gradient_template_file, 'Template');
Template = loaded_gradient.Template;

FrameShift = repmat(struct('sensor_id', NaN, 'base_xc_mm', NaN, ...
    'gradient_xc_mm', NaN, 'x_grid_shift_mm', NaN), numel(analysis_sensors), 1);

for is = 1:numel(analysis_sensors)
    sid = analysis_sensors(is);
    ib = find([loaded_base.Template.Sensor.sensor_id] == sid, 1);
    ig = find([Template.Sensor.sensor_id] == sid, 1);
    if isempty(ib) || isempty(ig)
        error('Missing CH%d in base or GradientXRange030 template.', sid);
    end

    base_xc = loaded_base.Template.Sensor(ib).xc;
    gradient_xc = Template.Sensor(ig).xc;
    x_shift = gradient_xc - base_xc;

    Template.Sensor(ig).x_grid = Template.Sensor(ig).x_grid + x_shift;
    if isfield(Template.Sensor(ig), 'x_domain') && ~isempty(Template.Sensor(ig).x_domain)
        Template.Sensor(ig).x_domain = Template.Sensor(ig).x_domain + x_shift;
    end
    Template.Sensor(ig).xc_gradient_original = gradient_xc;
    Template.Sensor(ig).xc_base_frame = base_xc;
    Template.Sensor(ig).xc = base_xc;
    Template.Sensor(ig).x_frame = 'base-template DynamicMap frame';

    FrameShift(is).sensor_id = sid;
    FrameShift(is).base_xc_mm = base_xc;
    FrameShift(is).gradient_xc_mm = gradient_xc;
    FrameShift(is).x_grid_shift_mm = x_shift;
end

Template.FrameMode = 'GradientXRange030 waveform in base-template x frame';
Template.SourceGradientTemplate = gradient_template_file;
Template.SourceBaseTemplateForFrame = base_template_file;
Template.FrameShift = struct2table(FrameShift);

save(base_frame_template_file, 'Template', '-v7.3');

fprintf('Saved BaseFrame GradientXRange030 template:\n  %s\n', base_frame_template_file);
disp(Template.FrameShift);
