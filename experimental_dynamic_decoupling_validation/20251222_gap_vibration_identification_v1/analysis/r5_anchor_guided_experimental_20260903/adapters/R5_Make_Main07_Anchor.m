function anchor = R5_Make_Main07_Anchor(Template, bladeId, sensorId)
%R5_MAKE_MAIN07_ANCHOR Build a frozen low-speed multi-lap anchor query.
rows = Template.SensorBlade;
idx = find([rows.blade_id] == bladeId & [rows.sensor_id] == sensorId, 1);
if isempty(idx)
    error('R5:AnchorMissing', 'Main07 template has no B%d / CH%d anchor.', bladeId, sensorId);
end
r = rows(idx);
x = double(r.x_grid(:)); y = R5_VoltageToTemplateMv(r.v_grid(:), r);
mask = isfinite(x) & isfinite(y) & logical(r.valid_grid_mask(:));
if isfield(r,'domain_effective_mask') && numel(r.domain_effective_mask)==numel(mask)
    mask = mask & logical(r.domain_effective_mask(:));
end
anchor = struct('x_mm',x,'y_mv',y,'supportMask',mask, ...
    'blade_id',bladeId,'sensor_id',sensorId, ...
    'template_lap_range',[r.stable_window_start_lap r.stable_window_end_lap], ...
    'source_file',Template.sourceFile,'voltage_contract','(V-baseline)*1000 mV');
end
