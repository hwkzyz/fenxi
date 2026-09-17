function R=run_inv_log_2_structured_main(lowInput,highMap,templateLib,cfg,mode)
%RUN_INV_LOG_2_STRUCTURED_MAIN Known-structure gap-vibration identification.
% Supported modes: single_sync, single_async, dual_sync_sync, dual_sync_async.
mainDir=fileparts(mfilename('fullpath'));localDir=fullfile(mainDir,'local_func');
if isempty(which('run_single_sync_voltage_vp')),addpath(localDir,'-begin');end
mode=lower(string(mode));allowed=["single_sync","single_async","dual_sync_sync","dual_sync_async"];
if ~any(mode==allowed),error('invlog2:UnknownStructuredMode','Unsupported structured mode: %s',mode);end
if mode=="dual_sync_sync"
    cfg.dualSyncSearchStrategy=get_field(cfg,'dualSyncSearchStrategy',"funnel_v1");
    cfg.funnelAllPairsFirstJointStep=get_field(cfg,'funnelAllPairsFirstJointStep',true);
    cfg.funnelSecondJointPairCount=get_field(cfg,'funnelSecondJointPairCount',24);
    cfg.funnelFullRefineCount=get_field(cfg,'funnelFullRefineCount',3);
end
nRev=numel(unique(highMap.rev_v(:)));minRev=get_field(cfg,'structuredMinRevolutions',8);
requireMin=get_field(cfg,'structuredRequireMinRevolutions',true);
if requireMin&&nRev<minRev
    error('invlog2:InsufficientRevolutions',...
        'Structured identification requires at least %d revolutions; received %d.',minRev,nRev);
end
% Direct single-/dual-frequency calls with a raw or aggregated low-speed
% record use the same frozen adaptive-SG calibration as the public entry.
if string(get_field(cfg,'route30ForwardModel',"low_increment"))=="low_increment" && ...
        ~(isstruct(lowInput)&&isfield(lowInput,'gHat')&&isfield(lowInput,'dx0')) && ...
        ~(isfield(templateLib,'fixedPathIncrement')&&templateLib.fixedPathIncrement)
    cfg.route30ForwardModel="low_increment";
    cfg.route30LowTemplateMethod="adaptive_sg";
    cfg.route30SupportAware=get_field(cfg,'route30SupportAware',true);
    cfg.route30FrequencyStructureMode=mode;
    C=calibrate_inv_log_2_low_speed(lowInput,templateLib,cfg);
    R=run_inv_log_2_high_with_calibration(C,highMap,cfg);
    R.structured_mode=mode;R.low_speed_state=C.lowState;
    R.forward_model_mode="low_template_increment";
    return;
end
state=build_low_state(lowInput,templateLib,cfg);state.calibration_source="low_speed_waveform";
switch mode
    case "single_sync"
        R=run_single_sync_voltage_vp(highMap,templateLib,cfg,state);
    case "single_async"
        R=run_single_async_voltage_vp(highMap,templateLib,cfg,state);
    case "dual_sync_sync"
        R=run_dual_sync_voltage_vp(highMap,templateLib,cfg,state);
    case "dual_sync_async"
        R=run_sync_async_voltage_vp(highMap,templateLib,cfg,state);
end
R.structured_mode=mode;R.low_speed_state=state;R.observed_revolutions=nRev;
R.minimum_required_revolutions=minRev;R.information_requirement_met=nRev>=minRev;
confidenceRev=get_field(cfg,'structuredConfidenceRevolutions',16);
maxRev=max(confidenceRev,get_field(cfg,'structuredMaximumRevolutions',32));
R.minimum_confidence_revolutions=confidenceRev;R.maximum_recommended_revolutions=maxRev;
isAmbiguous=isfield(R,'identification_status')&&R.identification_status=="ambiguous";
R.requires_more_observations=isAmbiguous&&nRev<maxRev;
R.identification_acceptable=isfield(R,'identification_confident')&&R.identification_confident;
if R.identification_acceptable
    R.observation_action="accept";
    R.reported_frequency=R.f_id;
elseif R.requires_more_observations
    R.observation_action="acquire_more_data";
    R.reported_frequency=nan(size(R.f_id));
else
    R.observation_action="report_nonunique";
    R.reported_frequency=nan(size(R.f_id));
end
R.provisional_frequency=R.f_id;
if isfield(R,'eo_id')
    if R.identification_acceptable,R.reported_order=R.eo_id;else,R.reported_order=nan(size(R.eo_id));end
elseif isfield(R,'sync_eo')
    if R.identification_acceptable,R.reported_order=R.sync_eo;else,R.reported_order=NaN;end
end
if isfield(templateLib,'fixedPathIncrement')&&templateLib.fixedPathIncrement
    R.forward_model_mode="low_template_increment";
else
    R.forward_model_mode="absolute";
end
end

function state=build_low_state(lowInput,lib,cfg)
if isstruct(lowInput)&&isfield(lowInput,'gHat')&&isfield(lowInput,'dx0')
    state=lowInput;return;
end
% A low-speed path-increment model already carries the calibrated absolute
% gap and spatial registration. Re-estimating them against that model from
% the same low-speed waveform is degenerate because the zero-increment path
% is independent of the query gap; use the calibration state directly.
if isfield(lib,'fixedPathIncrement') && lib.fixedPathIncrement && ...
        isfield(lib,'pathCal') && isstruct(lib.pathCal)
    state=struct('gHat',lib.pathCal.g0,'dx0',lib.pathCal.tau,...
        'bestCost',0,'calibration_source',"low_speed_path_calibration",...
        'low_speed_state',true);
    return;
end
if isstruct(lowInput)&&isfield(lowInput,'mapped')
    map=lowInput.mapped;
elseif isstruct(lowInput)&&isfield(lowInput,'x_v')&&isfield(lowInput,'V_a')
    map=lowInput;
else
    low=aggregate_low_speed_template(lowInput,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid);
    map=low.mapped;
end
state=estimate_highspeed_static_gap_raw(map,lib,cfg);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
