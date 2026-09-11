function Audit = R5_Audit_SensorWiseDg(resultFile)
%R5_AUDIT_SENSORWISEDG Audit independent per-sensor dynamic gap corrections.
% This is a read-only audit of a newly generated R5 Main06 result.  A result
% file must be supplied explicitly; historical files are never auto-selected.
if nargin < 1 || isempty(resultFile)
    error('R5:AuditInput','Pass the newly generated Main06 MAT result explicitly.');
end
S = load(resultFile,'Result');
assert(isfield(S,'Result'),'R5:AuditInput','MAT file lacks Result.');
R = S.Result;
assert(isfield(R,'MethodContract'),'R5:AuditContract','Missing MethodContract.');
assert(isfield(R.MethodContract,'sensor_wise_dynamic_gap') && R.MethodContract.sensor_wise_dynamic_gap, ...
    'R5:AuditContract','Result is not a sensor-wise dg result.');
assert(~R.MethodContract.pc1_dynamic_state, ...
    'R5:AuditContract','Dynamic PC1 state is not allowed in the formal R5 result.');
assert(isfield(R,'WindowResult') && ~isempty(R.WindowResult), ...
    'R5:AuditData','Result has no WindowResult records.');
nSensor = numel(R.WindowResult(1).modelFits.gap_only.deltaGapMm);
dgNames = arrayfun(@(k)sprintf('dg%d_mm',k),1:nSensor,'UniformOutput',false);
rows = repmat(struct('parameter','', 'n',0, 'mean_mm',NaN, 'std_mm',NaN, ...
    'min_mm',NaN, 'max_mm',NaN, 'bound_fraction',NaN, ...
    'min_boundary_margin_mm',NaN, 'median_bound_width_mm',NaN, ...
    'jump_p95_mm',NaN),0,1);
sensorIds = 1:nSensor;
if isfield(R,'cfg') && isfield(R.cfg,'analysisSensors')
    sensorIds = double(R.cfg.analysisSensors(:).');
end
for i = 1:numel(dgNames)
    xAll = arrayfun(@(w) double(w.modelFits.gap_only.deltaGapMm(i)), R.WindowResult(:));
    lower = nan(size(xAll)); upper = nan(size(xAll));
    for iw = 1:numel(R.WindowResult)
        W = R.WindowResult(iw);
        if ~isfield(W,'DeltaGapProjection') || ~isfield(W.DeltaGapProjection,'AllTable')
            continue;
        end
        P = W.DeltaGapProjection.AllTable;
        required = {'sensorId','EO','deltaGapLowerMm','deltaGapUpperMm'};
        if ~all(ismember(required,P.Properties.VariableNames)), continue; end
        selectedEO = double(W.modelFits.gap_only.EO);
        idx = find(double(P.sensorId)==sensorIds(i) & double(P.EO)==round(selectedEO),1,'first');
        if ~isempty(idx)
            lower(iw) = double(P.deltaGapLowerMm(idx));
            upper(iw) = double(P.deltaGapUpperMm(idx));
        end
    end
    valid = isfinite(xAll); x = xAll(valid);
    validBound = valid & isfinite(lower) & isfinite(upper) & upper > lower;
    assert(sum(validBound)==sum(valid),'R5:AuditBounds', ...
        'Missing actual local dg bounds for %s in one or more formal windows.',dgNames{i});
    margin = min(xAll(validBound)-lower(validBound),upper(validBound)-xAll(validBound));
    width = upper(validBound)-lower(validBound);
    nearBoundary = margin <= max(1e-8,0.02*width);
    d = abs(diff(xAll)); d = d(isfinite(d));
    q = struct('parameter',dgNames{i}, 'n',numel(x), 'mean_mm',NaN, ...
        'std_mm',NaN, 'min_mm',NaN, 'max_mm',NaN, ...
        'bound_fraction',NaN, 'min_boundary_margin_mm',NaN, ...
        'median_bound_width_mm',NaN, 'jump_p95_mm',NaN);
    if ~isempty(x)
        q.mean_mm = mean(x); q.std_mm = std(x); q.min_mm = min(x); q.max_mm = max(x);
        q.bound_fraction = mean(nearBoundary);
        q.min_boundary_margin_mm = min(margin);
        q.median_bound_width_mm = median(width);
        if ~isempty(d), q.jump_p95_mm = prctile(d,95); end
    end
    rows = [rows; q]; %#ok<AGROW>
end
Audit = struct('resultFile',resultFile,'status','pass','parameters',struct2table(rows), ...
    'method',R.method,'boundaryDefinition','actual_selected_EO_sensor_window_bounds');
typicalHalfWidth = 0.5*median(Audit.parameters.median_bound_width_mm,'omitnan');
if any(Audit.parameters.bound_fraction > 0.25) || ...
        any(Audit.parameters.jump_p95_mm > max(typicalHalfWidth,eps))
    Audit.status = 'scientific_flag';
end
disp(Audit.parameters);
end
