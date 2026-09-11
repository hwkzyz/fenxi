function window = R5_Make_Main08_Window(DynamicMap, windowId, sensorId, Template)
%R5_MAKE_MAIN08_WINDOW Convert one frozen Main08 sensor stream to R5 schema.
W = DynamicMap.Window;
idx = find([W.window_id] == windowId, 1);
if isempty(idx), error('R5:WindowMissing','No DynamicMap window %d.',windowId); end
S = W(idx).Sensor;
is = find([S.sensor_id] == sensorId,1);
if isempty(is), error('R5:SensorMissing','Window %d lacks CH%d.',windowId,sensorId); end
s = S(is);
rows = Template.SensorBlade;
it = find([rows.sensor_id] == sensorId, 1);
if isempty(it), error('R5:TemplateSensorMissing','Main07 template lacks CH%d.',sensorId); end
vMv = R5_VoltageToTemplateMv(s.V(:), rows(it));
window = struct('window_id',W(idx).window_id,'lap_range',W(idx).lap_range, ...
    'rpm_hz',W(idx).rot_freq_mean_hz,'rotFreqHz',W(idx).rot_freq_mean_hz, ...
    'sensor_id',sensorId, ...
    'x_mm',double(s.x_rel(:)),'t_s',double(s.t(:)),'v_mv',vMv, ...
    'validMask',isfinite(s.x_rel(:)) & isfinite(s.t(:)) & isfinite(vMv), ...
    'voltage_contract','(V-baseline)*1000 mV');
end
