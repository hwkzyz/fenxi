function C=calibrate_inv_log_2_low_speed(lowInput,templateLib,cfg)
%CALIBRATE_INV_LOG_2_LOW_SPEED Build one reusable low-speed calibration.
if string(get_field(cfg,'route30ForwardModel',"absolute"))~="low_increment"
    error('invlog2:CalibrationMode',...
        'Reusable low-speed calibration requires route30ForwardModel="low_increment".');
end
mainDir=fileparts(mfilename('fullpath'));localDir=fullfile(mainDir,'local_func');
if isempty(which('build_low_speed_templates_adaptive')),addpath(localDir,'-begin');end
if isstruct(lowInput)&&isfield(lowInput,'mapped')&&isfield(lowInput,'templateLow')
    low=lowInput;source="aggregated_low_speed_waveform";
elseif isstruct(lowInput)&&(isfield(lowInput,'x_v')||isfield(lowInput,'gHat'))
    error('invlog2:MissingLowSpeedTemplate',...
        'Mapped samples or a gap state alone cannot define a reusable forward template.');
else
    method=lower(string(get_field(cfg,'route30LowTemplateMethod',"adaptive_sg")));
    if method=="adaptive_sg"
        opts=get_field(cfg,'route30LowTemplateOptions',struct());
        if get_field(cfg,'route30SupportAware',false)&&~isfield(opts,'mapWindowRatio')
            if ~isfield(lowInput,'spatial_window_ratio')||isempty(lowInput.spatial_window_ratio)
                error('support:MissingAcquisitionWindow','Support-aware raw calibration requires declared low-speed spatial_window_ratio metadata.');
            end
            opts.mapWindowRatio=lowInput.spatial_window_ratio;
        end
        low=build_low_speed_templates_adaptive(lowInput,cfg.alpha_k,cfg.R_tip,...
            templateLib.domain,templateLib.xGrid,'cv_sg',opts);
    elseif method=="legacy"
        low=aggregate_low_speed_template(lowInput,cfg.alpha_k,cfg.R_tip,...
            templateLib.domain,templateLib.xGrid);
    else
        error('invlog2:UnknownLowTemplateMethod',...
            'Unsupported cfg.route30LowTemplateMethod: %s',method);
    end
    source="raw_low_speed_waveform";
end
state=estimate_highspeed_static_gap_raw(low.mapped,templateLib,cfg);
opts=get_field(cfg,'route30PathCalibrationOptions',struct());
if ~isfield(opts,'xDomain'),opts.xDomain=templateLib.domain;end
if ~isfield(opts,'fitMask'),opts.fitMask=[true false false false];end
if ~isfield(opts,'fitGainOffset'),opts.fitGainOffset=false;end
[pathCal,pathFit]=calibrate_low_speed_path(templateLib,low,opts);
model=build_path_template_model(templateLib,low.templateLow,pathCal);
state.g_low_hat=pathCal.g0;state.gHat=pathCal.g0;state.dx0=0;
state.g_parameter_type="delta_gap";state.calibration_source=source;
lowTemplate=low;if isfield(lowTemplate,'mapped'),lowTemplate=rmfield(lowTemplate,'mapped');end
C=struct('version',"adaptive_low_speed_calibration_v1",'templateModel',model,...
    'lowState',state,'lowTemplate',lowTemplate,'pathCal',pathCal,...
    'pathCalibrationFit',pathFit,'low_speed_input_type',source,...
    'low_speed_template_method',string(get_field(low,'method',"legacy_turn_interpolation")),...
    'forward_model_mode',"low_increment",...
    'measurement_geometry',struct('alpha_k',cfg.alpha_k,'R_tip',cfg.R_tip,...
    'domain',templateLib.domain),...
    'template_options',get_field(cfg,'route30LowTemplateOptions',struct()));
C.low_support_observation=summarize_low_support(low.mapped,model,cfg);
if get_field(cfg,'route30SupportAware',false)
    physical=get_field(cfg,'route30PhysicalDomain',struct());
    if isempty(fieldnames(physical))
        physical=default_physical_domain(model,cfg);
    end
    require_physical_contract(physical);
    C.version="support_aware_low_speed_calibration_v2";
    C.physical_domain=physical;
    C.support_aware=true;
else
    C.support_aware=false;
end
if isfield(low,'selected_window_mm')
    C.low_speed_selected_window_mm=low.selected_window_mm;
    C.low_speed_selected_span_bins=low.selected_span;
    C.low_speed_template_cv=low.cv_candidates;
end
end

function p=default_physical_domain(model,cfg)
g=model.pathCache.gGrid(:).';
dxBounds=get_field(cfg,'route30DefaultDxBoundsMm',[-.03 .03]);
aMax=get_field(cfg,'route30DefaultAmplitudeMaxMm',[.25 .25]);
p=struct('dx_bounds_mm',dxBounds,'amplitude_max_mm',aMax,...
    'gap_bounds_mm',sort([min(g) max(g)]),...
    'source',"default support contract derived after low-speed path calibration");
end
function S=summarize_low_support(map,model,cfg)
ids=unique(map.S_v(:));limits=nan(numel(ids),2);
for i=1:numel(ids),q=map.S_v(:)==ids(i);limits(i,:)=[min(map.x_v(q)) max(map.x_v(q))];end
margin=derive_support_numeric_margin(model,cfg);
S=struct('sensor_ids',ids,'per_sensor_observed_mm',limits,...
    'common_observed_mm',[max(limits(:,1)) min(limits(:,2))],...
    'map_window_ratio',get_field(map,'window_ratio',NaN),...
    'numeric_margin_audit',margin);
end
function require_physical_contract(p)
required={'dx_bounds_mm','amplitude_max_mm','gap_bounds_mm','source'};
for i=1:numel(required)
    if ~isfield(p,required{i})||isempty(p.(required{i}))
        error('support:MissingPhysicalContract','cfg.route30PhysicalDomain.%s is required.',required{i});
    end
end
if strlength(string(p.source))==0,error('support:MissingPhysicalSource','Physical-domain source must be documented.');end
end
function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
