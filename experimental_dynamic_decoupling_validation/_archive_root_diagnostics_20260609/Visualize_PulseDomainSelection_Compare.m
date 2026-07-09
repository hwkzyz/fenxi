%% Visualize pulse/domain selection strategies
% Compare the actual dynamic points used by:
%   1) all pulses + soft domain + soft margin 0.20 mm
%   2) single main pulse + hard domain + soft margin 0 mm
%
% The script is diagnostic only. It does not modify any Step03/Step07 output.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
figRoot = fullfile(thisDir, 'selection_diagnostics', 'pulse_domain_selection');
if exist(figRoot, 'dir') ~= 7
    mkdir(figRoot);
end
doExport = parse_logical_env_local('SELECTION_DIAG_EXPORT', false);

windowOverride = strtrim(getenv('SELECTION_DIAG_WINDOWS'));
if isempty(windowOverride)
    windowsToPlot = [1, 7, 8, 9];
else
    windowsToPlot = sscanf(windowOverride, '%d').';
end

cases = struct([]);
cases(end+1).name = '20250527_Main18';
cases(end).rotDir = fullfile(thisDir, '20250527_low_speed_rotating_calibration');
cases(end).sensorIds = [1, 3, 6];
cases(end).templateFile = fullfile(cases(end).rotDir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S136_GradientXRange030_20250527.mat');
cases(end).dynamicFile = fullfile(cases(end).rotDir, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S136_SlidingWindows_20250527.mat');

cases(end+1).name = '20251222_Main18';
cases(end).rotDir = fullfile(thisDir, '20251222_low_speed_rotating_calibration');
cases(end).sensorIds = [1, 2, 3];
cases(end).templateFile = fullfile(cases(end).rotDir, 'output', 'templates', ...
    'Template_LowSpeedRotating_B1_S123_GradientXRange030_20251222.mat');
cases(end).dynamicFile = fullfile(cases(end).rotDir, 'output', 'dynamic_maps', ...
    'DynamicMap_B1_S123_SlidingWindows_GradientXRange030_20251222.mat');

strategy(1).tag = 'all_soft_m020';
strategy(1).label = 'all + soft, margin=0.20';
strategy(1).pulseMode = 'all';
strategy(1).domainMode = 'soft';
strategy(1).softMarginMm = 0.20;

strategy(2).tag = 'single_adaptive_hard_m000';
strategy(2).label = 'single + adaptive hard';
strategy(2).pulseMode = 'single';
strategy(2).domainMode = 'hard';
strategy(2).softMarginMm = 0.0;

method = default_method_local();

summaryRows = [];
for ic = 1:numel(cases)
    C = cases(ic);
    if exist(C.templateFile, 'file') ~= 2
        error('Missing template file: %s', C.templateFile);
    end
    if exist(C.dynamicFile, 'file') ~= 2
        error('Missing dynamic map file: %s', C.dynamicFile);
    end
    St = load(C.templateFile, 'Template');
    Sd = load(C.dynamicFile, 'DynamicMap');
    Template = St.Template;
    DynamicMap = Sd.DynamicMap;

    figDir = fullfile(figRoot, C.name);
    if doExport && exist(figDir, 'dir') ~= 7
        mkdir(figDir);
    end

    validWindows = windowsToPlot(windowsToPlot >= 1 & windowsToPlot <= numel(DynamicMap.Window));
    for iw = validWindows
        Wmap = DynamicMap.Window(iw);
        [rowsX, rowsTime] = plot_window_selection_local(C, Template, Wmap, strategy, method, figDir, doExport);
        summaryRows = [summaryRows; rowsX(:); rowsTime(:)]; %#ok<AGROW>
    end
end

summaryFile = fullfile(figRoot, 'PulseDomainSelection_PointCounts.csv');
if doExport && ~isempty(summaryRows)
    writetable(struct2table(summaryRows), summaryFile);
    fprintf('Saved summary:\n  %s\n', summaryFile);
end
if doExport
    fprintf('Saved figures under:\n  %s\n', figRoot);
else
    disp('Figures are left open in MATLAB. Set SELECTION_DIAG_EXPORT=1 to save PNG/PDF/CSV.');
end

function method = default_method_local()
method.weightFloor = 0.05;
method.domainMarginMm = 0.02;
method.amplitudeLimitMm = 0.50;
method.dxLimitMm = 0.20;
method.coverageSafetyMarginMm = 0.05;
method.queryGuardMm = method.amplitudeLimitMm + method.dxLimitMm + method.coverageSafetyMarginMm;
method.queryGuardMode = 'adaptive';
method.queryGuardQuantile = 95;
method.queryGuardSafetyMm = method.coverageSafetyMarginMm;
method.queryGuardMinMm = 0.12;
method.queryGuardMaxMm = method.queryGuardMm;
method.defaultSensorThreshold = 0.5;
method.dynamicEffectiveMode = 'gradient';
method.dynamicTemplateGradientMinRatio = 0.08;
method.dynamicTimeGradientMinRatio = 0.15;
method.dynamicPeakQuantile = 85;
method.fontName = 'Times New Roman';
method.fontSize = 8.5;
method.labelFontSize = 9;
end

function [rowsX, rowsTime] = plot_window_selection_local(C, Template, Wmap, strategy, method, figDir, doExport)
nSensor = numel(C.sensorIds);
nStrategy = numel(strategy);
rowsX = repmat(empty_row_local(), nSensor * nStrategy, 1);
rowsTime = rowsX;

data = cell(nSensor, nStrategy);
for is = 1:nSensor
    sid = C.sensorIds(is);
    Tpl = get_template_sensor_local(Template, sid);
    D = get_dynamic_sensor_local(Wmap, sid);
    for im = 1:nStrategy
        data{is, im} = evaluate_selection_local(D, Tpl, strategy(im), method);
    end
end

figX = figure('Name', sprintf('%s W%02d x-domain selection', C.name, Wmap.window_id), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 24 6.5 + 3.2*nSensor]);
tl = tiledlayout(figX, nSensor, nStrategy, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s, window %d, lap [%g %g], x-domain view', ...
    C.name, Wmap.window_id, Wmap.lap_range(1), Wmap.lap_range(2)), ...
    'FontName', method.fontName, 'FontSize', method.labelFontSize, 'Interpreter', 'none');

rowIdx = 0;
for is = 1:nSensor
    sid = C.sensorIds(is);
    Tpl = get_template_sensor_local(Template, sid);
    D = get_dynamic_sensor_local(Wmap, sid);
    for im = 1:nStrategy
        rowIdx = rowIdx + 1;
        R = data{is, im};
        ax = nexttile(tl);
        hold(ax, 'on');
        plot(ax, Tpl.x_grid(:), (Tpl.v_grid(:) - Tpl.baseline) * 1000, ...
            'Color', [0.05 0.05 0.05], 'LineWidth', 1.0, 'DisplayName', 'template');
        scatter(ax, D.x_rel(:), (D.V(:) - Tpl.baseline) * 1000, 9, ...
            [0.78 0.78 0.78], 'filled', 'DisplayName', 'raw dynamic');
        scatter(ax, D.x_rel(R.finalMask), (D.V(R.finalMask) - Tpl.baseline) * 1000, ...
            15, R.weight(R.finalMask), 'filled', 'DisplayName', 'used points');
        xline(ax, Tpl.x_domain(1), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.75);
        xline(ax, Tpl.x_domain(2), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.75);
        if strcmpi(strategy(im).domainMode, 'hard')
            xline(ax, R.querySafeDomain(1), ':', 'Color', [0.85 0.1 0.1], 'LineWidth', 0.9);
            xline(ax, R.querySafeDomain(2), ':', 'Color', [0.85 0.1 0.1], 'LineWidth', 0.9);
        end
        colormap(ax, parula(128));
        clim(ax, [0 1]);
        box(ax, 'on');
        set(ax, 'TickDir', 'in', 'FontName', method.fontName, 'FontSize', method.fontSize);
        xlabel(ax, 'x relative to template center (mm)', 'FontName', method.fontName, ...
            'FontSize', method.labelFontSize);
        ylabel(ax, 'Voltage - baseline (mV)', 'FontName', method.fontName, ...
            'FontSize', method.labelFontSize);
        title(ax, sprintf('CH%d, %s, q=%.3f, N=%d/%d', sid, strategy(im).label, ...
            R.queryGuardMm, nnz(R.finalMask), numel(R.finalMask)), ...
            'FontName', method.fontName, 'FontSize', method.fontSize, 'Interpreter', 'none');
        if is == 1 && im == nStrategy
            cb = colorbar(ax);
            cb.Label.String = 'point weight';
            cb.FontName = method.fontName;
            cb.FontSize = method.fontSize;
        end
        rowsX(rowIdx) = make_row_local(C.name, Wmap.window_id, 'x_domain', sid, strategy(im), R);
    end
end
if doExport
    export_pair_local(figX, fullfile(figDir, sprintf('%s_W%02d_xdomain_selection', C.name, Wmap.window_id)));
    close(figX);
end

figT = figure('Name', sprintf('%s W%02d time selection', C.name, Wmap.window_id), ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 24 6.5 + 3.2*nSensor]);
tl = tiledlayout(figT, nSensor, nStrategy, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s, window %d, lap [%g %g], time view', ...
    C.name, Wmap.window_id, Wmap.lap_range(1), Wmap.lap_range(2)), ...
    'FontName', method.fontName, 'FontSize', method.labelFontSize, 'Interpreter', 'none');

rowIdx = 0;
for is = 1:nSensor
    sid = C.sensorIds(is);
    Tpl = get_template_sensor_local(Template, sid);
    D = get_dynamic_sensor_local(Wmap, sid);
    for im = 1:nStrategy
        rowIdx = rowIdx + 1;
        R = data{is, im};
        ax = nexttile(tl);
        hold(ax, 'on');
        tRel = D.t(:) - min(D.t(:));
        vMv = (D.V(:) - Tpl.baseline) * 1000;
        plot(ax, tRel, vMv, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
        scatter(ax, tRel(R.pulseMask), vMv(R.pulseMask), 9, [0.35 0.55 0.9], ...
            'filled', 'DisplayName', 'pulse mask');
        scatter(ax, tRel(R.finalMask), vMv(R.finalMask), 18, R.weight(R.finalMask), ...
            'filled', 'DisplayName', 'used points');
        yline(ax, R.thresholdMv, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.75);
        colormap(ax, parula(128));
        clim(ax, [0 1]);
        box(ax, 'on');
        set(ax, 'TickDir', 'in', 'FontName', method.fontName, 'FontSize', method.fontSize);
        xlabel(ax, 't in window (s)', 'FontName', method.fontName, 'FontSize', method.labelFontSize);
        ylabel(ax, 'Voltage - baseline (mV)', 'FontName', method.fontName, 'FontSize', method.labelFontSize);
        title(ax, sprintf('CH%d, %s, pulses=%d, N=%d', sid, strategy(im).label, ...
            R.pulseSegmentCount, nnz(R.finalMask)), ...
            'FontName', method.fontName, 'FontSize', method.fontSize, 'Interpreter', 'none');
        if is == 1 && im == nStrategy
            cb = colorbar(ax);
            cb.Label.String = 'point weight';
            cb.FontName = method.fontName;
            cb.FontSize = method.fontSize;
        end
        rowsTime(rowIdx) = make_row_local(C.name, Wmap.window_id, 'time', sid, strategy(im), R);
    end
end
if doExport
    export_pair_local(figT, fullfile(figDir, sprintf('%s_W%02d_time_selection', C.name, Wmap.window_id)));
    close(figT);
end
end

function R = evaluate_selection_local(D, Tpl, S, method)
xRaw = D.x_rel(:);
tRaw = D.t(:);
vRaw = D.V(:);
thetaRaw = D.theta(:);
threshold = resolve_threshold_local(Tpl, method);

f0 = interp1(Tpl.x_grid(:), Tpl.v_grid(:), xRaw, 'pchip', NaN);
fx0 = interp1(Tpl.x_grid(:), Tpl.dv_dx(:), xRaw, 'pchip', NaN);
xGridDomain = [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))];

if strcmpi(S.pulseMode, 'all')
    [maskPulse, pulseSegmentCount] = isolate_all_pulses_local(vRaw, threshold);
else
    [maskPulse, pulseSegmentCount] = isolate_main_pulse_local(vRaw, threshold);
end

maskEffective = build_dynamic_effective_mask_local(xRaw, vRaw, tRaw, Tpl, threshold, method);
maskDomain = xRaw >= Tpl.x_domain(1) + method.domainMarginMm & ...
             xRaw <= Tpl.x_domain(2) - method.domainMarginMm;
maskFinite = isfinite(f0) & isfinite(fx0) & isfinite(vRaw) & isfinite(thetaRaw);
mask = maskPulse & maskEffective & maskDomain & maskFinite;

queryGuardMm = resolve_query_guard_mm_local(Tpl, xRaw, vRaw, mask, method);
querySafeDomain = [xGridDomain(1) + queryGuardMm, xGridDomain(2) - queryGuardMm];
maskQuerySafe = xRaw >= querySafeDomain(1) & xRaw <= querySafeDomain(2);
if strcmpi(S.domainMode, 'hard')
    safeMask = mask & maskQuerySafe;
    if nnz(safeMask) >= 8
        mask = safeMask;
    end
end
if nnz(mask) < 8
    mask = maskEffective & maskDomain & maskFinite;
    if strcmpi(S.domainMode, 'hard')
        safeMask = mask & maskQuerySafe;
        if nnz(safeMask) >= 8
            mask = safeMask;
        end
    end
end

weight = compute_point_weight_local(xRaw, tRaw, vRaw, Tpl, S, method, queryGuardMm);
R = struct();
R.finalMask = mask;
R.pulseMask = maskPulse;
R.effectiveMask = maskEffective;
R.domainMask = maskDomain;
R.querySafeMask = maskQuerySafe;
R.weight = weight;
R.thresholdMv = (threshold - Tpl.baseline) * 1000;
R.pulseSegmentCount = pulseSegmentCount;
R.querySafeDomain = querySafeDomain;
R.queryGuardMm = queryGuardMm;
R.xSelectedMin = min(xRaw(mask), [], 'omitnan');
R.xSelectedMax = max(xRaw(mask), [], 'omitnan');
R.pointCount = nnz(mask);
R.rawCount = numel(mask);
R.pulseCount = nnz(maskPulse);
R.effectiveCount = nnz(maskEffective);
R.domainCount = nnz(maskDomain);
R.querySafeCount = nnz(maskQuerySafe);
end

function weight = compute_point_weight_local(x, t, v, Tpl, S, method, queryGuardMm)
wEdge = build_edge_weight_local(t, v, method.weightFloor);
if isfield(Tpl, 'bin_weight') && ~isempty(Tpl.bin_weight)
    wTemplate = interp1(Tpl.x_grid(:), Tpl.bin_weight(:), x(:), 'linear', method.weightFloor);
else
    wTemplate = ones(size(x));
end
wDomain = build_domain_soft_weight_local(x, [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))], ...
    S.softMarginMm, method.weightFloor);
wQuery = build_query_guard_soft_weight_local(x, [min(Tpl.x_grid(:)), max(Tpl.x_grid(:))], ...
    queryGuardMm, S.softMarginMm, method.weightFloor);
wGradient = build_template_gradient_weight_local(Tpl, x, S.domainMode, method.weightFloor);
weight = max(method.weightFloor, wEdge(:) .* wTemplate(:) .* wDomain(:) .* wQuery(:) .* wGradient(:));
if max(weight, [], 'omitnan') > 0
    weight = max(method.weightFloor, weight ./ max(weight, [], 'omitnan'));
end
weight(~isfinite(weight)) = method.weightFloor;
end

function threshold = resolve_threshold_local(Tpl, method)
if isfield(Tpl, 'threshold') && ~isempty(Tpl.threshold) && isfinite(Tpl.threshold)
    threshold = Tpl.threshold;
else
    threshold = Tpl.baseline + method.defaultSensorThreshold;
end
end

function queryGuardMm = resolve_query_guard_mm_local(Tpl, x, v, baseMask, method)
queryGuardMm = method.queryGuardMm;
if ~strcmpi(method.queryGuardMode, 'adaptive') || nnz(baseMask) < 8
    return;
end
xSel = x(baseMask);
vSel = v(baseMask);
xStat = invert_template_voltage_local(Tpl, vSel, xSel);
uApp = abs(xSel(:) - xStat(:));
uApp = uApp(isfinite(uApp));
if isempty(uApp)
    return;
end
q = min(max(method.queryGuardQuantile, 0), 100);
adaptiveGuard = prctile(uApp, q) + method.queryGuardSafetyMm;
queryGuardMm = min(max(adaptiveGuard, method.queryGuardMinMm), method.queryGuardMaxMm);
end

function xStat = invert_template_voltage_local(Tpl, v, xRef)
xGrid = Tpl.x_grid(:);
vGrid = Tpl.v_grid(:);
xStat = NaN(size(v));
for i = 1:numel(v)
    vv = v(i);
    diffV = vGrid - vv;
    crossingX = [];
    exactIdx = find(abs(diffV) <= 1e-10);
    if ~isempty(exactIdx)
        crossingX = xGrid(exactIdx);
    end
    for k = 1:numel(diffV)-1
        if ~isfinite(diffV(k)) || ~isfinite(diffV(k+1))
            continue;
        end
        if diffV(k) == 0 || diffV(k) * diffV(k+1) > 0
            continue;
        end
        denom = vGrid(k+1) - vGrid(k);
        if abs(denom) < eps
            continue;
        end
        alpha = (vv - vGrid(k)) / denom;
        crossingX(end+1, 1) = xGrid(k) + alpha * (xGrid(k+1) - xGrid(k)); %#ok<AGROW>
    end
    if isempty(crossingX)
        [~, idx] = min(abs(diffV));
        xStat(i) = xGrid(idx);
    else
        [~, idx] = min(abs(crossingX - xRef(i)));
        xStat(i) = crossingX(idx);
    end
end
end

function [mask, segmentCount] = isolate_all_pulses_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
for i = 1:numel(starts)
    mask(idx(starts(i)):idx(ends(i))) = true;
end
segmentCount = numel(starts);
end

function [mask, segmentCount] = isolate_main_pulse_local(v, threshold)
mask = false(size(v));
idx = find(v > threshold);
if isempty(idx)
    [~, imax] = max(v);
    span = max(5, round(numel(v) * 0.15));
    mask(max(1, imax-span):min(numel(v), imax+span)) = true;
    segmentCount = double(any(mask));
    return;
end
jumps = find(diff(idx) > 3);
starts = [1; jumps + 1];
ends = [jumps; numel(idx)];
if isscalar(starts)
    mask(idx(starts(1)):idx(ends(1))) = true;
    segmentCount = 1;
    return;
end

segmentGap = idx(starts(2:end)) - idx(ends(1:end-1));
largeGapThreshold = max(10, round(0.02 * numel(v)));
pulseBreaks = find(segmentGap > largeGapThreshold);
groupStarts = [1; pulseBreaks + 1];
groupEnds = [pulseBreaks; numel(starts)];
for ig = 1:numel(groupStarts)
    segIds = groupStarts(ig):groupEnds(ig);
    bestSeg = segIds(1);
    bestPeak = -inf;
    for iseg = segIds
        seg = idx(starts(iseg):ends(iseg));
        peakVal = max(v(seg));
        if peakVal > bestPeak
            bestPeak = peakVal;
            bestSeg = iseg;
        end
    end
    mask(idx(starts(bestSeg)):idx(ends(bestSeg))) = true;
end
segmentCount = numel(groupStarts);
end

function mask = build_dynamic_effective_mask_local(x, v, t, Tpl, threshold, method)
finite = isfinite(x) & isfinite(v) & isfinite(t);
if ~strcmpi(method.dynamicEffectiveMode, 'gradient')
    mask = finite;
    return;
end
gTpl = zeros(size(x));
if isfield(Tpl, 'dv_dx') && ~isempty(Tpl.dv_dx)
    gTpl = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
end
gTplMax = max(gTpl(finite), [], 'omitnan');
if ~isfinite(gTplMax) || gTplMax <= 0
    maskTpl = finite;
else
    maskTpl = gTpl >= method.dynamicTemplateGradientMinRatio * gTplMax;
end
gTime = zeros(size(v));
if nnz(finite) >= 5 && range(t(finite)) > 0
    gTime(finite) = abs(gradient(v(finite), t(finite)));
end
gTimeMax = max(gTime(finite), [], 'omitnan');
if ~isfinite(gTimeMax) || gTimeMax <= 0
    maskTime = false(size(v));
else
    maskTime = gTime >= method.dynamicTimeGradientMinRatio * gTimeMax;
end
peakLevel = prctile(v(finite), min(max(method.dynamicPeakQuantile, 0), 100));
maskPeak = v >= max(threshold, peakLevel);
mask = finite & (maskTpl | maskTime | maskPeak) & v >= 0.5 * threshold;
if nnz(mask) < 8
    mask = finite & v >= threshold;
end
if nnz(mask) < 8
    mask = finite;
end
end

function wEdge = build_edge_weight_local(t, v, floorW)
if numel(v) < 3 || range(t) <= 0
    wEdge = ones(size(v));
    return;
end
dv = abs(gradient(v(:), t(:)));
if max(dv) > 0
    wEdge = dv ./ max(dv);
else
    wEdge = ones(size(dv));
end
wEdge = max(floorW, wEdge);
end

function wDomain = build_domain_soft_weight_local(x, xDomain, marginMm, floorW)
if marginMm <= 0
    wDomain = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
ratio = min(max(distToEdge ./ marginMm, 0), 1);
wDomain = floorW + (1 - floorW) .* ratio;
wDomain(~isfinite(wDomain)) = floorW;
end

function wQuery = build_query_guard_soft_weight_local(x, xDomain, queryGuardMm, marginMm, floorW)
if queryGuardMm <= 0 || marginMm <= 0
    wQuery = ones(size(x));
    return;
end
distToEdge = min(x(:) - xDomain(1), xDomain(2) - x(:));
softStart = max(queryGuardMm - marginMm, 0);
ratio = min(max((distToEdge - softStart) ./ max(marginMm, eps), 0), 1);
wQuery = floorW + (1 - floorW) .* ratio;
wQuery(~isfinite(wQuery)) = floorW;
end

function wGradient = build_template_gradient_weight_local(Tpl, x, domainMode, floorW)
if ~strcmpi(domainMode, 'soft') || ~isfield(Tpl, 'dv_dx') || isempty(Tpl.dv_dx)
    wGradient = ones(size(x));
    return;
end
g = abs(interp1(Tpl.x_grid(:), Tpl.dv_dx(:), x(:), 'linear', 0));
gMax = max(g, [], 'omitnan');
if ~isfinite(gMax) || gMax <= 0
    wGradient = ones(size(x));
    return;
end
wGradient = floorW + (1 - floorW) .* g ./ gMax;
wGradient(~isfinite(wGradient)) = floorW;
end

function Tpl = get_template_sensor_local(Template, sid)
idx = find([Template.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Template missing CH%d.', sid);
end
Tpl = Template.Sensor(idx);
end

function D = get_dynamic_sensor_local(Wmap, sid)
idx = find([Wmap.Sensor.sensor_id] == sid, 1);
if isempty(idx)
    error('Dynamic window %d missing CH%d.', Wmap.window_id, sid);
end
D = Wmap.Sensor(idx);
end

function export_pair_local(fig, baseFile)
exportgraphics(fig, [baseFile '.png'], 'Resolution', 220);
exportgraphics(fig, [baseFile '.pdf'], 'ContentType', 'vector');
end

function tf = parse_logical_env_local(name, defaultValue)
raw = lower(strtrim(getenv(name)));
if isempty(raw)
    tf = defaultValue;
    return;
end
tf = any(strcmp(raw, {'1', 'true', 'yes', 'on'}));
end

function row = empty_row_local()
row = struct('dataset', "", 'window_id', NaN, 'view', "", 'sensor_id', NaN, ...
    'strategy', "", 'pulse_mode', "", 'domain_mode', "", 'soft_margin_mm', NaN, ...
    'raw_count', NaN, 'pulse_count', NaN, 'pulse_segments', NaN, ...
    'effective_count', NaN, 'domain_count', NaN, 'query_safe_count', NaN, ...
    'selected_count', NaN, 'x_selected_min_mm', NaN, 'x_selected_max_mm', NaN);
end

function row = make_row_local(dataset, windowId, viewName, sid, S, R)
row = empty_row_local();
row.dataset = string(dataset);
row.window_id = windowId;
row.view = string(viewName);
row.sensor_id = sid;
row.strategy = string(S.tag);
row.pulse_mode = string(S.pulseMode);
row.domain_mode = string(S.domainMode);
row.soft_margin_mm = S.softMarginMm;
row.raw_count = R.rawCount;
row.pulse_count = R.pulseCount;
row.pulse_segments = R.pulseSegmentCount;
row.effective_count = R.effectiveCount;
row.domain_count = R.domainCount;
row.query_safe_count = R.querySafeCount;
row.selected_count = R.pointCount;
row.x_selected_min_mm = R.xSelectedMin;
row.x_selected_max_mm = R.xSelectedMax;
end
