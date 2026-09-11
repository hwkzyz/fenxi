%% Independent preview-level ablation of OPR angle mappings for 20250527 V1.
% This diagnostic never overwrites the formal template, calibration library,
% or identification results.  It uses the stored formal waveform previews as
% a common observation set and rebuilds low-speed templates for three modes:
%   original_full_revolution, adjacent_uniform_60, unequal_physical_slot.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
packageDir = fileparts(fileparts(thisDir));
validationDir = fileparts(packageDir);
foundationDir = fullfile(validationDir, '20250527_btt_data_foundation');

lowSourceFile = fullfile(foundationDir, 'output', 'step04_low_speed_template', ...
    'LowSpeed_Template_SourceData_20250527.mat');
lowFeatureFile = fullfile(foundationDir, 'output', 'step01_low_speed_reference', ...
    'LowSpeed_Features_20250527.mat');
dynamicOprFile = fullfile(packageDir, 'inputs', 'prepared', 'foundation', ...
    'step02_dynamic_btt', '20250526_2500-3500_t400', 'jiluOPR.mat');
formalResultFile = fullfile(packageDir, 'results', 'fixed_gap', ...
    '20250526_2500-3500_t400', 'FixedGap_B1_S136_T001p5s.mat');
formalTemplateFile = fullfile(packageDir, 'inputs', 'prepared', 'foundation', ...
    'step04_low_speed_template', 'Template_OPRCenterStd_LowSpeed_AllBlades_S136_20250527.mat');

requiredFiles = {lowSourceFile, lowFeatureFile, dynamicOprFile, formalResultFile, formalTemplateFile};
for i = 1:numel(requiredFiles)
    assert(isfile(requiredFiles{i}), 'Missing diagnostic input: %s', requiredFiles{i});
end

P = struct();
P.dataset = '20250527';
P.sensorIds = [1 3 6];
P.targetBlade = 1;
P.eventsPerRevolution = 6;
P.rTipMm = 62;
P.gridDxMm = 0.02;
P.templateHalfWidthMm = 3.35;
P.minBinCount = 3;
P.smoothSpan = 9;
P.vpTopK = 3;
P.amplitudeLimitMm = 0.80;
P.dxLimitMm = 1.00;
P.frequencySearchHz = [300 1000];
P.modeNames = ["original_full_revolution", "adjacent_uniform_60", ...
    "unequal_physical_slot"];

fprintf('Loading low-speed OPR events...\n');
L = load(lowFeatureFile, 'case_data');
lowOprTimes = double(L.case_data.opr_times(:));
clear L;

fprintf('Loading low-speed point cloud...\n');
S = load(lowSourceFile, 'point_cloud');
pointCloud = S.point_cloud;
clear S;

H = load(dynamicOprFile, 'jiluOPR');
highOprTimes = double(H.jiluOPR(:, 1));
clear H;

F = load(formalResultFile, 'Result');
FormalResult = F.Result;
clear F;

%% 1. Calibrate unequal slot angles and low/high cyclic alignment.
[slotAngleDeg, slotMadDeg, perRevTable] = calibrate_slot_angles_local( ...
    lowOprTimes, P.eventsPerRevolution);
slotAnchorDeg = [0, cumsum(slotAngleDeg(1:end-1))];
[slotOffset, alignmentScore, alignmentMargin, alignmentTable] = ...
    align_high_to_low_local(lowOprTimes, highOprTimes, P.eventsPerRevolution);

assert(numel(slotAngleDeg) == P.eventsPerRevolution);
assert(abs(sum(slotAngleDeg) - 360) < 1e-9);
assert(isfinite(slotOffset));

SlotTable = table((1:P.eventsPerRevolution).', slotAnchorDeg(:), ...
    slotAngleDeg(:), slotMadDeg(:), ...
    'VariableNames', {'PhysicalSlot','AnchorAngleDeg','IntervalAngleDeg','MADDeg'});
writetable(SlotTable, fullfile(thisDir, 'OPR6_SlotAngles_20250527.csv'));
writetable(perRevTable, fullfile(thisDir, 'OPR6_SlotAngles_PerRevolution_20250527.csv'));
writetable(alignmentTable, fullfile(thisDir, 'OPR6_DynamicSlotAlignment_20250527.csv'));

SlotCalibration = struct();
SlotCalibration.schemaVersion = 'OPR6_PREVIEW_DIAGNOSTIC_V1';
SlotCalibration.dataset = P.dataset;
SlotCalibration.intervalAngleDeg = slotAngleDeg;
SlotCalibration.anchorAngleDeg = slotAnchorDeg;
SlotCalibration.slotIndexOffsetHighToLow = slotOffset;
SlotCalibration.alignmentScoreDeg = alignmentScore;
SlotCalibration.alignmentMarginDeg = alignmentMargin;
SlotCalibration.lowSourceFile = lowFeatureFile;
SlotCalibration.highSourceFile = dynamicOprFile;
save(fullfile(thisDir, 'OPR6_SlotCalibration_20250527.mat'), ...
    'SlotCalibration', 'SlotTable', 'alignmentTable');

fprintf('Slot angles: %s deg\n', mat2str(slotAngleDeg, 7));
fprintf('High-to-low slot offset: %d, score %.6f deg, margin %.6f deg\n', ...
    slotOffset, alignmentScore, alignmentMargin);

%% 2. Rebuild three low-speed template sets from the same point cloud.
TemplateByMode = struct();
for im = 1:numel(P.modeNames)
    modeName = P.modeNames(im);
    fprintf('Building template mode: %s\n', modeName);
    TemplateByMode.(char(modeName)) = build_templates_local( ...
        pointCloud, lowOprTimes, slotAngleDeg, slotAnchorDeg, 0, modeName, P);
    publish_compatible_template_local(formalTemplateFile, ...
        TemplateByMode.(char(modeName)), SlotCalibration, thisDir, P);
end

%% 3. Re-identify every stored formal waveform preview.
rows = repmat(empty_result_row_local(), 0, 1);
nWindow = numel(FormalResult.WindowResult);
for iw = 1:nWindow
    W = FormalResult.WindowResult(iw);
    B = choose_preview_bundle_local(W);
    if isempty(B) || ~isfield(B, 't') || isempty(B.t)
        warning('Window %d has no stored waveform preview.', iw);
        continue;
    end
    for im = 1:numel(P.modeNames)
        modeName = P.modeNames(im);
        Tmode = TemplateByMode.(char(modeName));
        [xMapped, thetaMapped] = remap_high_preview_local( ...
            B, highOprTimes, slotAngleDeg, slotAnchorDeg, slotOffset, ...
            modeName, Tmode, P);
        fit = identify_preview_local(B, xMapped, thetaMapped, Tmode, P);
        row = empty_result_row_local();
        row.WindowID = iw;
        row.TimeStartSec = W.time_window_s(1);
        row.TimeEndSec = W.time_window_s(2);
        row.MappingMode = modeName;
        row.EO = fit.EO;
        row.FrequencyHz = fit.frequencyHz;
        row.AmplitudeMm = fit.A;
        row.PhaseRad = fit.phi;
        row.DxMm = fit.dx;
        row.PlainVoltageRmseMv = 1000 * fit.rmse;
        row.PointCount = fit.pointCount;
        row.FormalEO = W.Result.EO_id;
        row.FormalAmplitudeMm = W.Result.A_id;
        row.FormalPlainRmseMv = 1000 * W.Result.plain_voltage_rmse;
        rows(end + 1, 1) = row; %#ok<SAGROW>
        fprintf('W%02d %-25s EO=%2d A=%.5f RMSE=%.3f mV\n', ...
            iw, modeName, fit.EO, fit.A, 1000 * fit.rmse);
    end
end

Trend = struct2table(rows);
writetable(Trend, fullfile(thisDir, 'OPR6_AngleMapping_PreviewTrend_20250527.csv'));
Summary = summarize_modes_local(Trend, P.modeNames);
writetable(Summary, fullfile(thisDir, 'OPR6_AngleMapping_PreviewSummary_20250527.csv'));

Diagnostic = struct();
Diagnostic.settings = P;
Diagnostic.slotCalibration = SlotCalibration;
Diagnostic.trend = Trend;
Diagnostic.summary = Summary;
Diagnostic.templateByMode = TemplateByMode;
Diagnostic.note = ['Preview-level closed-loop screening only. Formal raw-point ' ...
    'results and gap calibration are not modified.'];
save(fullfile(thisDir, 'OPR6_AngleMapping_PreviewAnalysis_20250527.mat'), ...
    'Diagnostic', '-v7.3');

disp(Summary);


function [slotAngleDeg, slotMadDeg, T] = calibrate_slot_angles_local(t, epr)
t = t(:);
n = numel(t) - epr;
idx = (1:n).';
revPeriod = t(idx + epr) - t(idx);
interval = t(idx + 1) - t(idx);
angle = 360 * interval ./ revPeriod;
slot = mod(idx - 1, epr) + 1;
valid = isfinite(angle) & angle > 0 & angle < 120;
slotAngleDeg = nan(1, epr);
slotMadDeg = nan(1, epr);
for s = 1:epr
    x = angle(valid & slot == s);
    slotAngleDeg(s) = median(x, 'omitnan');
    slotMadDeg(s) = median(abs(x - median(x, 'omitnan')), 'omitnan');
end
slotAngleDeg = 360 * slotAngleDeg / sum(slotAngleDeg);
revolution = floor((idx - 1) / epr) + 1;
T = table(idx, revolution, slot, angle, valid, ...
    'VariableNames', {'EventIndex','Revolution','EventSlot','IntervalAngleDeg','Valid'});
end


function [bestOffset, bestScore, margin, T] = align_high_to_low_local(tLow, tHigh, epr)
lowFingerprint = interval_fingerprint_local(tLow, epr);
highFingerprint = interval_fingerprint_local(tHigh, epr);
score = nan(epr, 1);
for offset = 0:(epr - 1)
    % Physical low slot assigned to high event k is mod(k-1+offset,epr)+1.
    observedByPhysicalSlot = circshift(highFingerprint, [0, offset]);
    score(offset + 1) = sqrt(mean((observedByPhysicalSlot - lowFingerprint).^2));
end
[bestScore, idx] = min(score);
bestOffset = idx - 1;
ss = sort(score);
margin = ss(2) - ss(1);
T = table((0:epr-1).', score, ...
    'VariableNames', {'SlotIndexOffsetHighToLow','FingerprintRmseDeg'});
end


function fp = interval_fingerprint_local(t, epr)
t = t(:);
n = numel(t) - epr;
idx = (1:n).';
a = 360 * (t(idx + 1) - t(idx)) ./ (t(idx + epr) - t(idx));
fp = nan(1, epr);
for s = 1:epr
    fp(s) = median(a(mod(idx - 1, epr) + 1 == s), 'omitnan');
end
fp = 360 * fp / sum(fp);
end


function Tmode = build_templates_local(pc, oprTimes, slotAngles, slotAnchors, ...
    slotOffset, modeName, P)
sensorTemplates = repmat(struct('sensorId', NaN, 'xGrid', [], 'vGrid', [], ...
    'dvDx', [], 'xCenterRawMm', NaN, 'xDomain', [], 'pointCount', 0), ...
    numel(P.sensorIds), 1);
for is = 1:numel(P.sensorIds)
    sid = P.sensorIds(is);
    idx = find([pc.sensor_id] == sid & [pc.blade_id] == P.targetBlade, 1);
    assert(~isempty(idx), 'Missing low-speed point cloud for CH%d B%d.', sid, P.targetBlade);
    q = pc(idx);
    switch char(modeName)
        case 'original_full_revolution'
            xRaw = double(q.x_mm(:));
        otherwise
            [angleLocalDeg, ~] = map_times_by_adjacent_slots_local( ...
                double(q.t_s(:)), oprTimes, slotAngles, slotAnchors, ...
                slotOffset, modeName, P.eventsPerRevolution);
            xRaw = P.rTipMm * pi / 180 * angleLocalDeg;
    end
    pulse = double(q.pulse_index(:));
    v = double(q.v(:));
    pulseCenters = pulse_centers_local(xRaw, v, pulse);
    xCenter = median(pulseCenters, 'omitnan');
    x = xRaw - xCenter;
    [xGrid, vGrid, dvDx] = aggregate_template_local(x, v, P);
    sensorTemplates(is).sensorId = sid;
    sensorTemplates(is).xGrid = xGrid;
    sensorTemplates(is).vGrid = vGrid;
    sensorTemplates(is).dvDx = dvDx;
    sensorTemplates(is).xCenterRawMm = xCenter;
    sensorTemplates(is).xDomain = [max(min(xGrid), -P.templateHalfWidthMm), ...
        min(max(xGrid), P.templateHalfWidthMm)];
    sensorTemplates(is).pointCount = numel(x);
end
Tmode = struct('mode', modeName, 'sensor', sensorTemplates, ...
    'slotAngleDeg', slotAngles, 'slotAnchorDeg', slotAnchors, ...
    'slotOffsetLow', slotOffset);
end


function publish_compatible_template_local(sourceFile, Tmode, slotCal, outDir, P)
S = load(sourceFile, 'Template');
Template = S.Template;
for is = 1:numel(P.sensorIds)
    sid = P.sensorIds(is);
    idx = find([Template.SensorBlade.sensor_id] == sid & ...
        [Template.SensorBlade.blade_id] == P.targetBlade, 1);
    assert(~isempty(idx), 'Compatible template source is missing CH%d B%d.', sid, P.targetBlade);
    src = Tmode.sensor(is);
    Template.SensorBlade(idx).x_grid = src.xGrid;
    Template.SensorBlade(idx).v_grid = src.vGrid;
    Template.SensorBlade(idx).v_grid_raw = src.vGrid;
    Template.SensorBlade(idx).v_grid_baseline_removed = src.vGrid - min(src.vGrid);
    Template.SensorBlade(idx).dv_dx = src.dvDx;
    Template.SensorBlade(idx).weight_grid = ones(size(src.xGrid));
    Template.SensorBlade(idx).count_grid = ones(size(src.xGrid));
    Template.SensorBlade(idx).valid_grid_mask = true(size(src.xGrid));
    Template.SensorBlade(idx).domain_effective_mask = ...
        src.xGrid >= src.xDomain(1) & src.xGrid <= src.xDomain(2);
    Template.SensorBlade(idx).domain_mask = Template.SensorBlade(idx).domain_effective_mask;
    Template.SensorBlade(idx).x_domain = src.xDomain;
    Template.SensorBlade(idx).x_coordinate_frame = 'OPRCenterStd_centered';
    Template.SensorBlade(idx).x_center_subtracted = true;
    Template.SensorBlade(idx).xc = src.xCenterRawMm;
    Template.SensorBlade(idx).xc_mm = src.xCenterRawMm;
    Template.SensorBlade(idx).xc_abs = src.xCenterRawMm;
    Template.SensorBlade(idx).x_grid_abs = src.xGrid + src.xCenterRawMm;
    Template.SensorBlade(idx).x_domain_abs = src.xDomain + src.xCenterRawMm;
end
Template.Metadata.angle_map_mode = char(Tmode.mode);
Template.Metadata.slot_angle_deg = Tmode.slotAngleDeg;
Template.Metadata.slot_anchor_deg = Tmode.slotAnchorDeg;
Template.Metadata.slot_index_offset = slotCal.slotIndexOffsetHighToLow;
Template.Metadata.slot_alignment_score = slotCal.alignmentScoreDeg;
Template.Metadata.slot_alignment_margin = slotCal.alignmentMarginDeg;
Template.Metadata.diagnostic_only = true;
metadata = Template.Metadata; %#ok<NASGU>
file = fullfile(outDir, sprintf('Template_Diagnostic_%s_20250527.mat', char(Tmode.mode)));
save(file, 'Template', 'metadata', '-v7.3');
end


function centers = pulse_centers_local(x, v, pulse)
ids = unique(pulse(isfinite(pulse)));
centers = nan(numel(ids), 1);
for i = 1:numel(ids)
    m = pulse == ids(i) & isfinite(x) & isfinite(v);
    xx = x(m); vv = v(m);
    if numel(xx) < 8
        continue;
    end
    sv = sort(vv);
    base = median(sv(1:max(3, floor(0.20 * numel(sv)))));
    amp = max(vv) - base;
    use = vv >= base + 0.85 * amp;
    if nnz(use) < 2
        [~, k] = max(vv);
        centers(i) = xx(k);
    else
        w = max(vv(use) - (base + 0.85 * amp), eps);
        centers(i) = sum(xx(use) .* w) / sum(w);
    end
end
end


function [xGrid, vGrid, dvDx] = aggregate_template_local(x, v, P)
valid = isfinite(x) & isfinite(v) & abs(x) <= 6;
x = x(valid); v = v(valid);
edges = (-6 - P.gridDxMm/2):P.gridDxMm:(6 + P.gridDxMm/2);
[~, ~, bin] = histcounts(x, edges);
nBin = numel(edges) - 1;
count = accumarray(bin(bin > 0), 1, [nBin, 1], @sum, 0);
vMed = accumarray(bin(bin > 0), v(bin > 0), [nBin, 1], @median, NaN);
xGrid = 0.5 * (edges(1:end-1) + edges(2:end));
good = count >= P.minBinCount & isfinite(vMed);
assert(nnz(good) > 30, 'Insufficient valid template bins.');
vGrid = interp1(xGrid(good), vMed(good), xGrid, 'pchip', NaN);
inside = xGrid >= min(xGrid(good)) & xGrid <= max(xGrid(good));
span = min(P.smoothSpan, nnz(inside));
if mod(span, 2) == 0
    span = span - 1;
end
vGrid(inside) = smoothdata(vGrid(inside), 'sgolay', max(3, span));
finiteGrid = isfinite(vGrid);
xGrid = xGrid(finiteGrid);
vGrid = vGrid(finiteGrid);
dvDx = gradient(vGrid, P.gridDxMm);
xGrid = xGrid(:); vGrid = vGrid(:); dvDx = dvDx(:);
end


function B = choose_preview_bundle_local(W)
B = [];
names = {'ExpandedBundlePreview','CoreBundlePreview','BundlePreview'};
for i = 1:numel(names)
    if isfield(W, names{i}) && isstruct(W.(names{i})) && ...
            isfield(W.(names{i}), 't') && ~isempty(W.(names{i}).t)
        B = W.(names{i});
        return;
    end
end
end


function [x, theta] = remap_high_preview_local(B, oprTimes, slotAngles, ...
    slotAnchors, slotOffset, modeName, Tmode, P)
switch char(modeName)
    case 'original_full_revolution'
        x = double(B.x(:));
        theta = double(B.theta(:));
        return;
end
[angleLocalDeg, theta] = map_times_by_adjacent_slots_local( ...
    double(B.t(:)), oprTimes, slotAngles, slotAnchors, slotOffset, ...
    modeName, P.eventsPerRevolution);
xRaw = P.rTipMm * pi / 180 * angleLocalDeg;
x = nan(size(xRaw));
sidVec = double(B.sensor_id(:));
for is = 1:numel(P.sensorIds)
    sid = P.sensorIds(is);
    m = sidVec == sid;
    x(m) = xRaw(m) - Tmode.sensor(is).xCenterRawMm;
end
end


function [angleLocalDeg, thetaRot] = map_times_by_adjacent_slots_local( ...
    t, oprTimes, slotAngles, slotAnchors, slotOffset, modeName, epr)
t = t(:); oprTimes = oprTimes(:);
prev = discretize(t, [-inf; oprTimes(2:end); inf]);
prev = max(1, min(prev, numel(oprTimes) - 1));
frac = (t - oprTimes(prev)) ./ (oprTimes(prev + 1) - oprTimes(prev));
eventInRev = mod(prev - 1, epr) + 1;
physicalSlot = mod(eventInRev - 1 + slotOffset, epr) + 1;
switch char(modeName)
    case 'adjacent_uniform_60'
        intervalDeg = 60 * ones(size(t));
        anchorDeg = 60 * (physicalSlot - 1);
    case 'unequal_physical_slot'
        intervalDeg = slotAngles(physicalSlot).';
        anchorDeg = slotAnchors(physicalSlot).';
    otherwise
        error('Unsupported adjacent-slot mode: %s', modeName);
end
angleLocalDeg = frac .* intervalDeg;
revIndex = floor((prev - 1) / epr);
thetaRot = 2*pi*revIndex + anchorDeg*pi/180 + frac .* intervalDeg*pi/180;
invalid = t < oprTimes(1) | t > oprTimes(end) | ~isfinite(frac);
angleLocalDeg(invalid) = NaN;
thetaRot(invalid) = NaN;
end


function fit = identify_preview_local(B, x, theta, Tmode, P)
v = double(B.v(:));
sid = double(B.sensor_id(:));
mask = isfinite(x) & isfinite(theta) & isfinite(v) & ismember(sid, P.sensorIds);
if isfield(B, 'main_pulse_mask')
    mask = mask & logical(B.main_pulse_mask(:));
end
for is = 1:numel(P.sensorIds)
    m = sid == P.sensorIds(is);
    d = Tmode.sensor(is).xDomain;
    mask(m) = mask(m) & x(m) >= d(1) + P.amplitudeLimitMm + 0.05 & ...
        x(m) <= d(2) - P.amplitudeLimitMm - 0.05;
end
x = x(mask); theta = theta(mask); v = v(mask); sid = sid(mask);
assert(numel(v) > 100, 'Too few preview points after common query-safe mask.');

rotFreqHz = double(B.rot_freq_mean_hz);
eoMin = max(1, ceil(P.frequencySearchHz(1) / rotFreqHz));
eoMax = floor(P.frequencySearchHz(2) / rotFreqHz);
eoList = eoMin:eoMax;
seed = repmat(struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx', NaN, 'rmse', inf), ...
    numel(eoList), 1);

[t0, g0] = evaluate_template_stack_local(Tmode, sid, x);
y = v - t0;
for i = 1:numel(eoList)
    eo = eoList(i);
    X = [-g0, -g0 .* sin(eo * theta), -g0 .* cos(eo * theta)];
    beta = X \ y;
    dx = beta(1);
    A = hypot(beta(2), beta(3));
    phi = atan2(beta(3), beta(2));
    [rmse, ~] = objective_local([A, phi, dx], eo, x, theta, v, sid, Tmode, P);
    seed(i) = struct('EO', eo, 'A', A, 'phi', phi, 'dx', dx, 'rmse', rmse);
end
[~, order] = sort([seed.rmse]);
order = order(1:min(P.vpTopK, numel(order)));
best = struct('EO', NaN, 'A', NaN, 'phi', NaN, 'dx', NaN, 'rmse', inf);
opt = optimset('Display', 'off', 'MaxIter', 350, 'MaxFunEvals', 1400, ...
    'TolX', 1e-9, 'TolFun', 1e-10);
for k = order(:).'
    eo = seed(k).EO;
    p0 = [min(max(seed(k).A, 0), P.amplitudeLimitMm), seed(k).phi, ...
        min(max(seed(k).dx, -P.dxLimitMm), P.dxLimitMm)];
    fun = @(p) objective_local(p, eo, x, theta, v, sid, Tmode, P);
    p = fminsearch(fun, p0, opt);
    [rmse, pb] = objective_local(p, eo, x, theta, v, sid, Tmode, P);
    if rmse < best.rmse
        best = struct('EO', eo, 'A', pb(1), 'phi', pb(2), ...
            'dx', pb(3), 'rmse', rmse);
    end
end
fit = best;
fit.frequencyHz = fit.EO * rotFreqHz;
fit.pointCount = numel(v);
end


function [rmse, p] = objective_local(pRaw, eo, x, theta, v, sid, Tmode, P)
p = pRaw(:).';
p(1) = min(max(abs(p(1)), 0), P.amplitudeLimitMm);
p(2) = mod(p(2) + pi, 2*pi) - pi;
p(3) = min(max(p(3), -P.dxLimitMm), P.dxLimitMm);
xq = x - p(3) - p(1) * sin(eo * theta + p(2));
pred = evaluate_template_stack_local(Tmode, sid, xq);
r = v - pred;
valid = isfinite(r);
if nnz(valid) < 100
    rmse = 1e3;
else
    rmse = sqrt(mean(r(valid).^2));
end
outside = nnz(~valid) / numel(valid);
rmse = rmse + 10 * outside;
end


function [v, g] = evaluate_template_stack_local(Tmode, sid, x)
v = nan(size(x));
g = nan(size(x));
for is = 1:numel(Tmode.sensor)
    m = sid == Tmode.sensor(is).sensorId;
    v(m) = interp1(Tmode.sensor(is).xGrid, Tmode.sensor(is).vGrid, x(m), 'pchip', NaN);
    if nargout > 1
        g(m) = interp1(Tmode.sensor(is).xGrid, Tmode.sensor(is).dvDx, x(m), 'linear', NaN);
    end
end
end


function row = empty_result_row_local()
row = struct('WindowID', NaN, 'TimeStartSec', NaN, 'TimeEndSec', NaN, ...
    'MappingMode', "", 'EO', NaN, 'FrequencyHz', NaN, 'AmplitudeMm', NaN, ...
    'PhaseRad', NaN, 'DxMm', NaN, 'PlainVoltageRmseMv', NaN, ...
    'PointCount', NaN, 'FormalEO', NaN, 'FormalAmplitudeMm', NaN, ...
    'FormalPlainRmseMv', NaN);
end


function S = summarize_modes_local(T, modes)
rows = repmat(struct('MappingMode', "", 'WindowCount', NaN, ...
    'ModeEO', NaN, 'ModeEOCount', NaN, 'MeanFrequencyHz', NaN, ...
    'MeanAmplitudeMm', NaN, 'MeanPlainVoltageRmseMv', NaN, ...
    'MedianPlainVoltageRmseMv', NaN), numel(modes), 1);
for i = 1:numel(modes)
    Q = T(T.MappingMode == modes(i), :);
    rows(i).MappingMode = modes(i);
    rows(i).WindowCount = height(Q);
    if ~isempty(Q)
        rows(i).ModeEO = mode(Q.EO);
        rows(i).ModeEOCount = nnz(Q.EO == rows(i).ModeEO);
        rows(i).MeanFrequencyHz = mean(Q.FrequencyHz, 'omitnan');
        rows(i).MeanAmplitudeMm = mean(Q.AmplitudeMm, 'omitnan');
        rows(i).MeanPlainVoltageRmseMv = mean(Q.PlainVoltageRmseMv, 'omitnan');
        rows(i).MedianPlainVoltageRmseMv = median(Q.PlainVoltageRmseMv, 'omitnan');
    end
end
S = struct2table(rows);
end
