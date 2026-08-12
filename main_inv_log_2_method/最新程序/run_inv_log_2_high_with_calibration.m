function R=run_inv_log_2_high_with_calibration(C,highMap,cfg)
%RUN_INV_LOG_2_HIGH_WITH_CALIBRATION Reuse one low-speed calibration.
required={'templateModel','lowState','pathCal','pathCalibrationFit'};
for i=1:numel(required)
    if ~isfield(C,required{i}),error('invlog2:InvalidCalibration',...
            'Calibration is missing field %s.',required{i});end
end
supportContract=[];
if isfield(C,'support_aware')&&C.support_aware
    if ~isfield(C,'physical_domain'),error('support:MissingPhysicalContract','Support-aware calibration has no physical domain.');end
    if ~isfield(C,'low_support_observation'),error('support:MissingLowObservation','Support-aware calibration has no observed low-speed support.');end
    lowProxy=struct('x_v',reshape(C.low_support_observation.per_sensor_observed_mm.',[],1),...
        'S_v',repelem(C.low_support_observation.sensor_ids(:),2));
    supportContract=derive_support_aware_domain(highMap,lowProxy,C.templateModel,...
        C.physical_domain,cfg);
    if ~supportContract.gap_domain_supported
        error('support:GapDomainOutsideLibrary','Physical gap domain is outside the cached static-response support.');
    end
    [highMap,supportContract]=restrict_map_to_support_contract(highMap,supportContract);
    cfg.supportContract=supportContract;
    cfg.route30PhysicalDomain=C.physical_domain;
end
R=run_inv_log_2_main_method(highMap,C.templateModel,cfg,C.lowState);
R.dg_used=R.g_used-C.pathCal.g0;R.g_low_used=C.pathCal.g0;
R.gap_parameter_type="delta_gap";R.pathCal=C.pathCal;
R.pathCalibrationFit=C.pathCalibrationFit;R.low_speed_state=C.lowState;
R.low_speed_input_type=C.low_speed_input_type;R.forward_model_mode="low_increment";
R.low_speed_template_method=C.low_speed_template_method;
if ~isempty(supportContract),R.support_contract=supportContract;end
if isfield(C,'low_speed_selected_window_mm')
    R.low_speed_selected_window_mm=C.low_speed_selected_window_mm;
    R.low_speed_selected_span_bins=C.low_speed_selected_span_bins;
    R.low_speed_template_cv=C.low_speed_template_cv;
end
end
