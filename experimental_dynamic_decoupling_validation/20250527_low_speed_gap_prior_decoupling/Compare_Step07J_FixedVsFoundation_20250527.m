%% Compare Step07J gap-aware trends with the no-gap foundation Step05 trend.
% This diagnostic does not run Step07J by default. It reads existing trend
% CSV files, compares the synchronized no-gap branch against foundation
% Step05, and also reports how gap_only/gap_tilt move away from that
% baseline. Set STEP07J_FVF_STEP07J_TREND_FILE or
% STEP07J_FVF_FOUNDATION_TREND_FILE to override either input. Set
% STEP07J_FVF_ALLOW_MAIN_FALLBACK=1 to compare only the main gap_tilt trend
% when a method_compare trend has not been generated yet.

clear; clc; close all;

thisDir = fileparts(mfilename('fullpath'));
oldDir = pwd;
cleanupDir = onCleanup(@() cd(oldDir));
cd(thisDir);
addpath(thisDir);

C0 = CaseConfig();
rootDir = fileparts(thisDir);
outDir = fullfile(thisDir, 'outputs');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

[foundationDefault, foundationSource] = default_foundation_file_local(rootDir, C0);
foundationFile = env_or_default_local('STEP07J_FVF_FOUNDATION_TREND_FILE', foundationDefault);
if exist(foundationFile, 'file') ~= 2
    foundationFile = find_foundation_fallback_local(rootDir, C0);
end
if exist(foundationFile, 'file') ~= 2
    error('Missing foundation trend CSV. Set STEP07J_FVF_FOUNDATION_TREND_FILE or run foundation Step05 first.');
end

allowMainFallback = parse_logical_env_local('STEP07J_FVF_ALLOW_MAIN_FALLBACK', false);
step07jFile = resolve_step07j_trend_file_local(outDir, C0, allowMainFallback);
if exist(step07jFile, 'file') ~= 2
    error(['Missing Step07J trend CSV. Run the gap-folder Step07J comparison/main result first, ' ...
        'or set STEP07J_FVF_STEP07J_TREND_FILE.']);
end

FoundationRaw = readtable(foundationFile, 'TextType', 'string');
Step07JRaw = readtable(step07jFile, 'TextType', 'string');

foundationSpec = struct( ...
    'eo', 'EO_id', ...
    'freq', 'fn_id', ...
    'amp', 'A_id', ...
    'dx', 'dx_c_id', ...
    'rmse', 'weighted_voltage_rmse', ...
    'rmseScale', 1000, ...
    'rotFreq', 'rot_freq_mean_hz');
foundation = extract_series_local(FoundationRaw, 'foundation_step05', foundationSpec, foundationFile);

methodSpecs = step07j_method_specs_local();
methodSeries = cell(numel(methodSpecs), 1);
nMethodSeries = 0;
for i = 1:numel(methodSpecs)
    if is_method_available_local(Step07JRaw, methodSpecs(i))
        candidateSeries = extract_series_local( ...
            Step07JRaw, methodSpecs(i).method, methodSpecs(i), step07jFile);
        if series_has_data_local(candidateSeries)
            nMethodSeries = nMethodSeries + 1;
            methodSeries{nMethodSeries} = candidateSeries;
        end
    end
end
methodSeries = methodSeries(1:nMethodSeries);
if isempty(methodSeries)
    error('The Step07J trend CSV does not contain fixed/gap/tilt method columns: %s', step07jFile);
end
methodNames = cellfun(@(s) string(s.method), methodSeries);
if ~any(methodNames == "fixed") && ~allowMainFallback
    error(['The Step07J trend does not contain the fixed branch needed for fixed-vs-foundation. ' ...
        'Run Step09_Compare_Step07J_Methods_%s first, or set STEP07J_FVF_ALLOW_MAIN_FALLBACK=1.'], C0.dataset);
end

pairTables = cell(numel(methodSeries), 1);
for i = 1:numel(methodSeries)
    pairTables{i} = pair_with_foundation_local(foundation, methodSeries{i});
end
WindowCompare = vertcat(pairTables{:});
Summary = build_summary_local(foundation, methodSeries, WindowCompare);

suffix = sanitize_suffix_local(strtrim(getenv('STEP07J_FVF_SUFFIX')));
if isempty(suffix)
    if allowMainFallback && ~any(methodNames == "fixed")
        suffix = '_main_gaptilt_fallback_vs_foundation';
    else
        suffix = '_fixed_vs_foundation';
    end
end
summaryFile = fullfile(outDir, sprintf( ...
    'Compare_Step07J_FixedVsFoundation_Summary_%s_%s%s.csv', ...
    C0.dataset, C0.caseTag, suffix));
windowFile = fullfile(outDir, sprintf( ...
    'Compare_Step07J_FixedVsFoundation_Window_%s_%s%s.csv', ...
    C0.dataset, C0.caseTag, suffix));
figFile = fullfile(outDir, sprintf( ...
    'Compare_Step07J_FixedVsFoundation_Trend_%s_%s%s.png', ...
    C0.dataset, C0.caseTag, suffix));

writetable(Summary, summaryFile);
writetable(WindowCompare, windowFile);
plot_comparison_local(foundation, methodSeries, figFile, C0);

fprintf('\n=== Step07J fixed/gap vs foundation: %s %s ===\n', C0.dataset, C0.caseTag);
fprintf('Foundation: %s [%s]\n', foundationFile, foundationSource);
fprintf('Step07J:    %s\n', step07jFile);
disp(Summary);
fprintf('Saved summary:\n  %s\n', summaryFile);
fprintf('Saved window comparison:\n  %s\n', windowFile);
fprintf('Saved trend figure:\n  %s\n', figFile);

delete(cleanupDir);
cd(oldDir);

function [file, source] = default_foundation_file_local(rootDir, C)
source = 'dataset_default';
switch char(C.dataset)
    case '20241106'
        file = fullfile(rootDir, '20241106_btt_data_foundation', 'output', ...
            'step05_single_sync_direct_template_fullbundle_for_gap_step07j', '3000_3150', ...
            sprintf('Trend_Step05_FoundationMainPulseAdaptiveFixedJointEta_%s_%s.csv', ...
            C.caseTag, C.dataset));
    case '20250527'
        file = fullfile(rootDir, '20250527_btt_data_foundation', 'output', ...
            'step05_single_sync_direct_template_fullbundle_gap_step07j', '20250526_2500-3500_t400', ...
            sprintf('Trend_Step05_SingleSyncDirectTemplate_%s_%s.csv', C.caseTag, C.dataset));
    case '20251222'
        file = fullfile(rootDir, '20251222_btt_data_foundation', 'output', ...
            'step05_single_sync_direct_template_fullbundle_gap_step07j', '1000_2500_3500', ...
            sprintf('Trend_Step05_FixedEtaGradDispVP_%s_%s.csv', ...
            C.caseTag, C.dataset));
    otherwise
        error('Unsupported dataset: %s', C.dataset);
end
end

function value = env_or_default_local(name, defaultValue)
value = strtrim(getenv(name));
if isempty(value)
    value = defaultValue;
end
end

function file = find_foundation_fallback_local(rootDir, C)
foundationRoot = fullfile(rootDir, sprintf('%s_btt_data_foundation', C.dataset), 'output');
patterns = {
    sprintf('Trend_Step05_FoundationMainPulseAdaptiveFixedJointEta_%s_%s.csv', C.caseTag, C.dataset)
    sprintf('Trend_Step05_SingleSyncDirectTemplate_%s_%s.csv', C.caseTag, C.dataset)
    sprintf('Trend_Step05_*FixedJoint*%s*%s.csv', C.caseTag, C.dataset)
    };
file = '';
for i = 1:numel(patterns)
    hits = dir(fullfile(foundationRoot, '**', patterns{i}));
    if ~isempty(hits)
        file = pick_latest_file_local(hits);
        return;
    end
end
end

function file = resolve_step07j_trend_file_local(outDir, C, allowMainFallback)
override = strtrim(getenv('STEP07J_FVF_STEP07J_TREND_FILE'));
if ~isempty(override)
    file = override;
    return;
end
patterns = {
    sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s_method_compare.csv', C.dataset, C.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s_*method_compare.csv', C.dataset, C.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s_*nested_compare*.csv', C.dataset, C.caseTag)
    sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s_*compare*.csv', C.dataset, C.caseTag)
    };
if allowMainFallback
    patterns(end+1:end+2) = {
        sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s_main_gaptilt.csv', C.dataset, C.caseTag)
        sprintf('Step07J_NestedStaticWarp_VPFullWave_Trend_%s_%s_*main_gaptilt.csv', C.dataset, C.caseTag)
        };
end
file = '';
for i = 1:numel(patterns)
    hits = dir(fullfile(outDir, patterns{i}));
    if ~isempty(hits)
        file = pick_latest_file_local(hits);
        return;
    end
end
end

function file = pick_latest_file_local(hits)
[~, idx] = max([hits.datenum]);
file = fullfile(hits(idx).folder, hits(idx).name);
end

function specs = step07j_method_specs_local()
specs = struct( ...
    'method', {'direct_low_template', 'fixed', 'gap_only', 'gap_tilt'}, ...
    'eo', {'direct_EO', 'fixed_EO', 'gap_EO', 'tilt_EO'}, ...
    'freq', {'direct_frequency_hz', 'fixed_frequency_hz', 'gap_frequency_hz', 'tilt_frequency_hz'}, ...
    'amp', {'direct_amplitude_mm', 'fixed_amplitude_mm', 'gap_amplitude_mm', 'tilt_amplitude_mm'}, ...
    'dx', {'', 'fixed_dx_mm', 'gap_dx_mm', 'tilt_dx_mm'}, ...
    'rmse', {'direct_rmse_mV', 'fixed_rmse_mV', 'gap_rmse_mV', 'tilt_rmse_mV'}, ...
    'rmseScale', {1, 1, 1, 1}, ...
    'rotFreq', {'rot_freq_hz', 'rot_freq_hz', 'rot_freq_hz', 'rot_freq_hz'});
end

function tf = is_method_available_local(T, spec)
vars = T.Properties.VariableNames;
tf = ismember(spec.eo, vars) && ismember(spec.amp, vars) && ismember(spec.rmse, vars);
end

function tf = series_has_data_local(s)
tf = any(isfinite(s.EO)) || any(isfinite(s.amplitudeMm)) || any(isfinite(s.rmseMv));
end

function s = extract_series_local(T, methodName, spec, sourceFile)
n = height(T);
s = struct();
s.method = char(methodName);
s.sourceFile = sourceFile;
s.windowId = column_by_candidates_local(T, {'window_id'}, n);
s.lapStart = column_by_candidates_local(T, {'lap_start'}, n);
s.lapEnd = column_by_candidates_local(T, {'lap_end'}, n);
s.timeS = column_by_candidates_local(T, {'window_center_time_s'}, n);
s.rotFreqHz = column_by_candidates_local(T, {spec.rotFreq, 'rot_freq_hz', 'rot_freq_mean_hz'}, n);
s.EO = column_by_candidates_local(T, {spec.eo}, n);
s.frequencyHz = column_by_candidates_local(T, {spec.freq}, n);
s.amplitudeMm = column_by_candidates_local(T, {spec.amp}, n);
if isempty(spec.dx)
    s.dxMm = nan(n, 1);
else
    s.dxMm = column_by_candidates_local(T, {spec.dx}, n);
end
s.rmseMv = spec.rmseScale .* column_by_candidates_local(T, {spec.rmse}, n);
end

function x = column_by_candidates_local(T, candidates, n)
vars = T.Properties.VariableNames;
name = '';
for i = 1:numel(candidates)
    if ~isempty(candidates{i}) && ismember(candidates{i}, vars)
        name = candidates{i};
        break;
    end
end
if isempty(name)
    x = nan(n, 1);
else
    x = to_numeric_column_local(T.(name), n);
end
end

function x = to_numeric_column_local(raw, n)
if isnumeric(raw) || islogical(raw)
    x = double(raw(:));
elseif iscell(raw)
    x = str2double(string(raw(:)));
elseif isstring(raw)
    x = str2double(raw(:));
elseif iscategorical(raw)
    x = str2double(string(raw(:)));
else
    x = nan(n, 1);
end
if numel(x) ~= n
    y = nan(n, 1);
    y(1:min(n, numel(x))) = x(1:min(n, numel(x)));
    x = y;
end
end

function Pair = pair_with_foundation_local(foundation, series)
[fi, si, mode] = match_indices_local(foundation.windowId, series.windowId);
n = numel(fi);
method = repmat(string(series.method), n, 1);
matchMode = repmat(string(mode), n, 1);
pairIndex = (1:n).';
foundationWindowId = foundation.windowId(fi);
methodWindowId = series.windowId(si);
foundationLapStart = foundation.lapStart(fi);
foundationLapEnd = foundation.lapEnd(fi);
methodLapStart = series.lapStart(si);
methodLapEnd = series.lapEnd(si);
foundationEO = foundation.EO(fi);
methodEO = series.EO(si);
foundationFrequencyHz = foundation.frequencyHz(fi);
methodFrequencyHz = series.frequencyHz(si);
foundationAmplitudeMm = foundation.amplitudeMm(fi);
methodAmplitudeMm = series.amplitudeMm(si);
foundationDxMm = foundation.dxMm(fi);
methodDxMm = series.dxMm(si);
foundationRmseMv = foundation.rmseMv(fi);
methodRmseMv = series.rmseMv(si);
deltaFrequencyHz = methodFrequencyHz - foundationFrequencyHz;
deltaAmplitudeMm = methodAmplitudeMm - foundationAmplitudeMm;
deltaDxMm = methodDxMm - foundationDxMm;
deltaRmseMv = methodRmseMv - foundationRmseMv;
eoMatch = isfinite(foundationEO) & isfinite(methodEO) & foundationEO == methodEO;

Pair = table(method, matchMode, pairIndex, foundationWindowId, methodWindowId, ...
    foundationLapStart, foundationLapEnd, methodLapStart, methodLapEnd, ...
    foundationEO, methodEO, eoMatch, ...
    foundationFrequencyHz, methodFrequencyHz, deltaFrequencyHz, ...
    foundationAmplitudeMm, methodAmplitudeMm, deltaAmplitudeMm, ...
    foundationDxMm, methodDxMm, deltaDxMm, ...
    foundationRmseMv, methodRmseMv, deltaRmseMv);
end

function [fi, si, mode] = match_indices_local(fWindow, sWindow)
fFinite = isfinite(fWindow);
sFinite = isfinite(sWindow);
if any(fFinite) && any(sFinite)
    [common, fi, si] = intersect(fWindow(fFinite), sWindow(sFinite), 'stable');
    fMap = find(fFinite);
    sMap = find(sFinite);
    fi = fMap(fi);
    si = sMap(si);
    if numel(common) >= max(1, floor(0.5 * min(nnz(fFinite), nnz(sFinite))))
        mode = 'window_id';
        return;
    end
end
n = min(numel(fWindow), numel(sWindow));
fi = (1:n).';
si = (1:n).';
mode = 'row_index';
end

function Summary = build_summary_local(foundation, methodSeries, WindowCompare)
n = 1 + numel(methodSeries);
method = strings(n, 1);
sourceFile = strings(n, 1);
windowCount = zeros(n, 1);
eoDistribution = strings(n, 1);
dominantEO = nan(n, 1);
meanFrequencyHz = nan(n, 1);
stdFrequencyHz = nan(n, 1);
meanAmplitudeMm = nan(n, 1);
meanRmseMv = nan(n, 1);
medianRmseMv = nan(n, 1);
pairedWindowCount = nan(n, 1);
eoMatchFractionVsFoundation = nan(n, 1);
meanAbsAmplitudeDeltaVsFoundationMm = nan(n, 1);
medianAbsAmplitudeDeltaVsFoundationMm = nan(n, 1);
meanRmseDeltaVsFoundationMv = nan(n, 1);

allSeries = [{foundation}; methodSeries(:)];
for i = 1:n
    s = allSeries{i};
    method(i) = string(s.method);
    sourceFile(i) = string(s.sourceFile);
    finiteEO = isfinite(s.EO);
    windowCount(i) = nnz(finiteEO);
    eoDistribution(i) = eo_distribution_local(s.EO);
    dominantEO(i) = mode_finite_local(s.EO);
    meanFrequencyHz(i) = mean(s.frequencyHz, 'omitnan');
    stdFrequencyHz(i) = std(s.frequencyHz, 'omitnan');
    meanAmplitudeMm(i) = mean(s.amplitudeMm, 'omitnan');
    meanRmseMv(i) = mean(s.rmseMv, 'omitnan');
    medianRmseMv(i) = median(s.rmseMv, 'omitnan');
    if i > 1
        mask = WindowCompare.method == string(s.method);
        pairedWindowCount(i) = nnz(mask);
        eoMatchFractionVsFoundation(i) = mean(double(WindowCompare.eoMatch(mask)), 'omitnan');
        meanAbsAmplitudeDeltaVsFoundationMm(i) = mean(abs(WindowCompare.deltaAmplitudeMm(mask)), 'omitnan');
        medianAbsAmplitudeDeltaVsFoundationMm(i) = median(abs(WindowCompare.deltaAmplitudeMm(mask)), 'omitnan');
        meanRmseDeltaVsFoundationMv(i) = mean(WindowCompare.deltaRmseMv(mask), 'omitnan');
    end
end

Summary = table(method, sourceFile, windowCount, eoDistribution, dominantEO, ...
    meanFrequencyHz, stdFrequencyHz, meanAmplitudeMm, meanRmseMv, medianRmseMv, ...
    pairedWindowCount, eoMatchFractionVsFoundation, ...
    meanAbsAmplitudeDeltaVsFoundationMm, medianAbsAmplitudeDeltaVsFoundationMm, ...
    meanRmseDeltaVsFoundationMv);
end

function text = eo_distribution_local(v)
v = v(isfinite(v));
if isempty(v)
    text = "";
    return;
end
vals = unique(v(:)).';
parts = strings(1, numel(vals));
for i = 1:numel(vals)
    parts(i) = sprintf('EO%d=%d', vals(i), sum(v == vals(i)));
end
text = strjoin(parts, ', ');
end

function y = mode_finite_local(v)
v = v(isfinite(v));
if isempty(v)
    y = NaN;
    return;
end
vals = unique(v(:));
counts = zeros(size(vals));
for i = 1:numel(vals)
    counts(i) = sum(v == vals(i));
end
[~, idx] = max(counts);
y = vals(idx);
end

function plot_comparison_local(foundation, methodSeries, figFile, C)
fig = figure('Name', sprintf('%s %s Step07J vs foundation', C.dataset, C.caseTag), ...
    'Color', 'w', 'Visible', 'off', 'Units', 'centimeters', 'Position', [2, 2, 18, 16]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on; grid on; box on;
plot_series_local(foundation, foundation.EO, 'k-', 'foundation');
for i = 1:numel(methodSeries)
    s = methodSeries{i};
    plot_series_local(s, s.EO, style_for_method_local(s.method), s.method);
end
ylabel('EO'); title(sprintf('%s %s EO', C.dataset, C.caseTag)); legend('Location', 'best');

nexttile; hold on; grid on; box on;
plot_series_local(foundation, foundation.amplitudeMm, 'k-', 'foundation');
for i = 1:numel(methodSeries)
    s = methodSeries{i};
    plot_series_local(s, s.amplitudeMm, style_for_method_local(s.method), s.method);
end
ylabel('Amplitude (mm)'); title('Identified amplitude');

nexttile; hold on; grid on; box on;
plot_series_local(foundation, foundation.rmseMv, 'k-', 'foundation');
for i = 1:numel(methodSeries)
    s = methodSeries{i};
    plot_series_local(s, s.rmseMv, style_for_method_local(s.method), s.method);
end
xlabel('Window'); ylabel('Weighted RMSE (mV)'); title('Full waveform objective');

try
    exportgraphics(fig, figFile, 'Resolution', 300);
catch
    saveas(fig, figFile);
end
close(fig);
end

function plot_series_local(s, y, style, labelText)
x = s.windowId;
if all(~isfinite(x))
    x = (1:numel(y)).';
end
mask = isfinite(x) & isfinite(y);
if any(mask)
    plot(x(mask), y(mask), style, 'LineWidth', 1.2, 'MarkerSize', 3.8, ...
        'DisplayName', char(labelText));
end
end

function style = style_for_method_local(methodName)
switch char(methodName)
    case 'direct_low_template'
        style = '-.';
    case 'fixed'
        style = '-o';
    case 'gap_only'
        style = '-s';
    case 'gap_tilt'
        style = '-^';
    otherwise
        style = '--';
end
end

function tf = parse_logical_env_local(name, defaultValue)
txt = lower(strtrim(getenv(name)));
if isempty(txt)
    tf = defaultValue;
elseif any(strcmp(txt, {'1','true','yes','on'}))
    tf = true;
elseif any(strcmp(txt, {'0','false','no','off'}))
    tf = false;
else
    error('%s must be true/false.', name);
end
end

function suffix = sanitize_suffix_local(suffix)
if isempty(suffix)
    return;
end
if suffix(1) ~= '_'
    suffix = ['_', suffix];
end
suffix = regexprep(suffix, '[^\w\-]', '_');
end
