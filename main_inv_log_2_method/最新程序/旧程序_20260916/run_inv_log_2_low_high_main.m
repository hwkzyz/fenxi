function R = run_inv_log_2_low_high_main(lowSpeedInput, highMap, templateLib, cfg)
%RUN_INV_LOG_2_LOW_HIGH_MAIN End-to-end low/high-speed main entry.
%
% lowSpeedInput may be:
%   1) a raw low-speed waveform accepted by map_highspeed_to_space;
%   2) an aggregate_low_speed_template result containing .mapped;
%   3) an already mapped structure containing .x_v and .V_a;
%   4) an estimated low-speed state containing .gHat and .dx0.

if nargin < 4
    error('lowSpeedInput, highMap, templateLib, and cfg are required.');
end

mainDir = fileparts(mfilename('fullpath'));
localDir = fullfile(mainDir, 'local_func');
if isempty(which('estimate_highspeed_static_gap_raw'))
    addpath(localDir, '-begin');
end

% Structured single/dual routes use the frozen support-aware adaptive-SG
% calibration.  Keep the historical absolute route as the default for the
% legacy general router unless the caller explicitly selects low_increment.
if isfield(cfg,'route30ForwardModel') && ~isempty(cfg.route30ForwardModel)
    forwardMode=string(cfg.route30ForwardModel);
else
    structureMode0=lower(string(get_field(cfg,'route30FrequencyStructureMode',"general")));
    if ismember(structureMode0,["single_sync","single_async","dual_sync_sync","dual_sync_async"])
        forwardMode="low_increment";
    else
        forwardMode="absolute";
    end
end
if forwardMode == "low_increment"
    cfg.route30ForwardModel = "low_increment";
    cfg.route30LowTemplateMethod = "adaptive_sg";
    structureMode=lower(string(get_field(cfg,'route30FrequencyStructureMode',"general")));
    defaultSupportAware=ismember(structureMode,...
        ["single_sync","single_async","dual_sync_sync","dual_sync_async"]);
    cfg.route30SupportAware = get_field(cfg,'route30SupportAware',defaultSupportAware);
    C=calibrate_inv_log_2_low_speed(lowSpeedInput,templateLib,cfg);
    R=run_inv_log_2_high_with_calibration(C,highMap,cfg);
    return;
end
[lowState, source, low] = build_low_speed_state(lowSpeedInput, templateLib, cfg);
lowState.calibration_source = source;
R = run_inv_log_2_main_method(highMap, templateLib, cfg, lowState);
R.gap_parameter_type = "absolute_gap";
R.low_speed_state = lowState;
R.low_speed_input_type = source;
R.forward_model_mode = forwardMode;
if ~isempty(low)
    R.low_speed_template_method = string(get_field(low, 'method', "legacy_turn_interpolation"));
    if isfield(low, 'selected_window_mm')
        R.low_speed_selected_window_mm = low.selected_window_mm;
        R.low_speed_selected_span_bins = low.selected_span;
        R.low_speed_template_cv = low.cv_candidates;
    end
end
end

function [state, source, low] = build_low_speed_state(lowInput, templateLib, cfg)
low = [];
if isstruct(lowInput) && isfield(lowInput, 'gHat') && isfield(lowInput, 'dx0')
    state = lowInput;
    source = "provided_low_speed_state";
    return;
end

if isstruct(lowInput) && isfield(lowInput, 'mapped')
    low = lowInput;
    mapped = lowInput.mapped;
    source = "aggregated_low_speed_waveform";
elseif isstruct(lowInput) && isfield(lowInput, 'x_v') && isfield(lowInput, 'V_a')
    mapped = lowInput;
    source = "mapped_low_speed_waveform";
else
    required = {'alpha_k','R_tip'};
    for i = 1:numel(required)
        if ~isfield(cfg, required{i})
            error('cfg.%s is required to map a raw low-speed waveform.', required{i});
        end
    end
    if ~isfield(templateLib, 'xGrid')
        error('templateLib.xGrid is required to aggregate a raw low-speed waveform.');
    end
    low = aggregate_low_speed_template(lowInput, cfg.alpha_k, cfg.R_tip, ...
        templateLib.domain, templateLib.xGrid);
    mapped = low.mapped;
    source = "raw_low_speed_waveform";
end

state = estimate_highspeed_static_gap_raw(mapped, templateLib, cfg);
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
