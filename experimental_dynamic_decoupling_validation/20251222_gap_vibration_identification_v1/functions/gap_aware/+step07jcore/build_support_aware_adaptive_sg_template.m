function [Template, Diagnostics] = build_support_aware_adaptive_sg_template( ...
        pointCloud, Template, cfg)
%BUILD_SUPPORT_AWARE_ADAPTIVE_SG_TEMPLATE Grouped-lap adaptive SG template.
% Existing x coordinates and support/query-safe domains are frozen. Only
% the voltage template and its derivative are rebuilt from the same raw
% low-speed point cloud. PCHIP remains the downstream interpolation rule.

if nargin < 3 || isempty(cfg), cfg = struct(); end
candidateWindowMm = get_cfg_local(cfg, 'candidateWindowMm', ...
    [0.10 0.18 0.34 0.50 0.82 1.22 1.62 2.02 2.42]);
foldCount = get_cfg_local(cfg, 'foldCount', 5);
minBinCount = get_cfg_local(cfg, 'minBinCount', 5);
assert(isfield(Template, 'SensorBlade') && ~isempty(Template.SensorBlade), ...
    'Template.SensorBlade is required.');

entries = Template.SensorBlade;
bladeIds = unique([entries.blade_id], 'stable');
selectionRows = repmat(empty_selection_local(), numel(bladeIds), 1);
sensorRows = repmat(empty_sensor_local(), numel(entries), 1);
sensorRowCount = 0;

for ib = 1:numel(bladeIds)
    bladeId = bladeIds(ib);
    entryIndex = find([entries.blade_id] == bladeId);
    data = repmat(struct('x', [], 'v', [], 'lap', [], 'grid', [], ...
        'entryIndex', NaN, 'sensorId', NaN), numel(entryIndex), 1);
    dx = nan(numel(entryIndex), 1);
    for is = 1:numel(entryIndex)
        ie = entryIndex(is);
        sid = entries(ie).sensor_id;
        ip = find([pointCloud.sensor_id] == sid & ...
            [pointCloud.blade_id] == bladeId, 1);
        assert(~isempty(ip), 'Missing low-speed point cloud S%d B%d.', sid, bladeId);
        pc = pointCloud(ip);
        grid = entries(ie).x_grid(:);
        x = pc.x_mm(:) - entries(ie).xc;
        v = pc.v(:);
        lap = round(pc.lap_index(:));
        keep = isfinite(x) & isfinite(v) & isfinite(lap) & ...
            x >= min(grid) & x <= max(grid);
        if isfield(pc, 'stable_window_start_lap') && ...
                isfinite(pc.stable_window_start_lap) && ...
                isfield(pc, 'stable_window_end_lap') && ...
                isfinite(pc.stable_window_end_lap)
            keep = keep & lap >= pc.stable_window_start_lap & ...
                lap <= pc.stable_window_end_lap;
        end
        data(is).x = x(keep);
        data(is).v = v(keep);
        data(is).lap = lap(keep);
        data(is).grid = grid;
        data(is).entryIndex = ie;
        data(is).sensorId = sid;
        dx(is) = median(diff(grid));
    end

    allLaps = unique(vertcat(data.lap));
    actualFoldCount = min(foldCount, numel(allLaps));
    assert(actualFoldCount >= 2, 'At least two complete low-speed laps are required.');
    foldId = mod((1:numel(allLaps)) - 1, actualFoldCount) + 1;
    loss = nan(actualFoldCount, numel(data), numel(candidateWindowMm));
    for f = 1:actualFoldCount
        testLaps = allLaps(foldId == f);
        for is = 1:numel(data)
            train = ~ismember(data(is).lap, testLaps);
            test = ~train;
            if nnz(train) < 100 || nnz(test) < 20, continue; end
            raw = bin_template_local(data(is).x(train), data(is).v(train), ...
                data(is).grid, minBinCount);
            for ic = 1:numel(candidateWindowMm)
                span = physical_span_local(candidateWindowMm(ic), ...
                    median(diff(data(is).grid)), numel(data(is).grid));
                fit = sg_smooth_local(raw, span);
                loss(f, is, ic) = prediction_rmse_mv_local(data(is), test, fit);
            end
        end
    end
    z = reshape(loss, [], numel(candidateWindowMm));
    cvMean = mean(z, 1, 'omitnan');
    cvSe = std(z, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(z), 1));
    [~, choice] = min(cvMean);
    eligible = cvMean <= cvMean(choice) + 2 .* cvSe(choice);
    screenChoice = find(eligible, 1, 'last');
    assert(~isempty(choice) && isfinite(cvMean(choice)), ...
        'Grouped-lap CV failed for blade %d.', bladeId);

    selectionRows(ib).bladeId = bladeId;
    selectionRows(ib).sensorIds = join(string([data.sensorId]), ',');
    selectionRows(ib).foldCount = actualFoldCount;
    selectionRows(ib).selectedWindowMm = candidateWindowMm(choice);
    selectionRows(ib).selectedSpanBins = round(median(arrayfun(@(q) ...
        physical_span_local(candidateWindowMm(choice), q, inf), dx)));
    selectionRows(ib).cvRmseMv = cvMean(choice);
    selectionRows(ib).cvSeMv = cvSe(choice);
    selectionRows(ib).screenWindowMm = candidateWindowMm(screenChoice);

    for is = 1:numel(data)
        ie = data(is).entryIndex;
        raw = bin_template_local(data(is).x, data(is).v, data(is).grid, minBinCount);
        span = physical_span_local(candidateWindowMm(choice), ...
            median(diff(data(is).grid)), numel(data(is).grid));
        smooth = sg_smooth_local(raw, span);
        derivative = gradient(smooth, data(is).grid);
        entries(ie).v_grid = smooth;
        entries(ie).dv_dx = derivative;
        if isfield(entries, 'v_grid_baseline_removed')
            entries(ie).v_grid_baseline_removed = smooth - entries(ie).baseline;
        end
        if isfield(entries, 'amplitude')
            entries(ie).amplitude = max(smooth, [], 'omitnan') - entries(ie).baseline;
        end
        sensorRowCount = sensorRowCount + 1;
        sensorRows(sensorRowCount).sensorId = data(is).sensorId;
        sensorRows(sensorRowCount).bladeId = bladeId;
        sensorRows(sensorRowCount).pointCount = numel(data(is).x);
        sensorRows(sensorRowCount).lapCount = numel(unique(data(is).lap));
        sensorRows(sensorRowCount).selectedWindowMm = candidateWindowMm(choice);
        sensorRows(sensorRowCount).selectedSpanBins = span;
        sensorRows(sensorRowCount).supportLeftMm = min(data(is).grid);
        sensorRows(sensorRowCount).supportRightMm = max(data(is).grid);
        if isfield(entries(ie), 'x_query_safe_domain')
            sensorRows(sensorRowCount).querySafeLeftMm = ...
                entries(ie).x_query_safe_domain(1);
            sensorRows(sensorRowCount).querySafeRightMm = ...
                entries(ie).x_query_safe_domain(2);
        end
    end

    selectionRows(ib).candidateWindowMm = {candidateWindowMm(:).'};
    selectionRows(ib).candidateCvRmseMv = {cvMean(:).'};
    selectionRows(ib).candidateCvSeMv = {cvSe(:).'};
end

Template.SensorBlade = entries;
if isfield(Template, 'SchemaVersion') && strlength(string(Template.SchemaVersion)) > 0
    Template.SchemaVersion = [char(string(Template.SchemaVersion)), '/adaptiveSG-v1'];
else
    Template.SchemaVersion = 'adaptiveSG-v1';
end
Template.LowSpeed_Smoothing_Method = 'support_aware_grouped_lap_cv_sg_pchip';
Template.LowSpeed_Interpolation_Method = 'pchip_no_extrapolation';
Template.AdaptiveSG_Selection = struct2table(selectionRows);
if isfield(Template, 'Template_Settings')
    Template.Template_Settings.smoothing_method = ...
        'support_aware_grouped_lap_cv_sg';
    Template.Template_Settings.interpolation_method = 'pchip';
    Template.Template_Settings.sg_candidate_window_mm = candidateWindowMm;
    Template.Template_Settings.sg_grouped_fold_count = foldCount;
end
Diagnostics = struct();
Diagnostics.method = 'support_aware_grouped_lap_cv_sg_pchip';
Diagnostics.selection = struct2table(selectionRows);
Diagnostics.sensor = struct2table(sensorRows(1:sensorRowCount));
Diagnostics.candidateWindowMm = candidateWindowMm;
Diagnostics.cvGrouping = 'complete_low_speed_laps';
Diagnostics.selectionRule = 'minimum_grouped_lap_voltage_prediction_rmse';
Diagnostics.forwardInterpolation = 'pchip_no_extrapolation';
end


function y = bin_template_local(x, v, grid, minCount)
dx = median(diff(grid));
edges = [grid(1) - dx ./ 2; grid(1:end-1) + diff(grid) ./ 2; grid(end) + dx ./ 2];
bin = discretize(x, edges);
ok = isfinite(bin) & isfinite(v);
y = accumarray(bin(ok), v(ok), [numel(grid), 1], @median, NaN);
count = accumarray(bin(ok), 1, [numel(grid), 1], @sum, 0);
y(count < minCount) = NaN;
valid = isfinite(y);
assert(nnz(valid) >= 5, 'Insufficient valid low-speed bins.');
y = fillmissing(y, 'linear', 'SamplePoints', grid);
y = fillmissing(y, 'nearest');
end


function y = sg_smooth_local(y, span)
if span >= 5
    y = smoothdata(y(:), 'sgolay', span, 'Degree', min(3, span - 1));
else
    y = y(:);
end
end


function span = physical_span_local(widthMm, dxMm, n)
span = max(5, 2 .* floor((widthMm ./ dxMm) ./ 2) + 1);
if isfinite(n)
    span = min(span, n);
    if mod(span, 2) == 0, span = span - 1; end
end
end


function rmse = prediction_rmse_mv_local(data, test, fit)
vTest = data.v(test);
pred = interp1(data.grid, fit, data.x(test), 'pchip', NaN);
ok = isfinite(vTest) & isfinite(pred);
if nnz(ok) < 20
    rmse = NaN;
else
    rmse = 1000 .* sqrt(mean((vTest(ok) - pred(ok)).^2));
end
end


function value = get_cfg_local(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end


function row = empty_selection_local()
row = struct('bladeId', NaN, 'sensorIds', "", 'foldCount', NaN, ...
    'selectedWindowMm', NaN, 'selectedSpanBins', NaN, ...
    'cvRmseMv', NaN, 'cvSeMv', NaN, 'screenWindowMm', NaN, ...
    'candidateWindowMm', {{}}, 'candidateCvRmseMv', {{}}, ...
    'candidateCvSeMv', {{}});
end


function row = empty_sensor_local()
row = struct('sensorId', NaN, 'bladeId', NaN, 'pointCount', NaN, ...
    'lapCount', NaN, 'selectedWindowMm', NaN, 'selectedSpanBins', NaN, ...
    'supportLeftMm', NaN, 'supportRightMm', NaN, ...
    'querySafeLeftMm', NaN, 'querySafeRightMm', NaN);
end
