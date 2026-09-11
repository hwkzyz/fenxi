function Comparison = Compare_ParameterStudy_20241106(summaryFile, groupName)
%COMPARE_PARAMETERSTUDY_20241106 Compare completed isolated study cases.

studyDir = fileparts(mfilename('fullpath'));
addpath(studyDir);
Study = Config_ParameterStudy_20241106();
resultDir = Study.outputRoot;
if nargin < 1 || isempty(summaryFile)
    summaryFile = fullfile(resultDir, 'ParameterStudy_Summary_Latest.csv');
end
assert(isfile(summaryFile), 'Parameter-study summary does not exist: %s', summaryFile);
T = readtable(summaryFile, 'TextType', 'string');
T = T(strcmpi(T.status, "pass"), :);
if nargin >= 2 && ~isempty(groupName)
    T = T(strcmpi(T.group, string(groupName)), :);
end
assert(~isempty(T), 'No completed parameter-study cases match the selection.');

referenceIndex = find_reference_local(T);
referenceAmplitude = T.mean_amplitude_mm(referenceIndex);
referenceFrequency = T.mean_frequency_hz(referenceIndex);
referenceRmse = T.mean_plain_rmse_mV(referenceIndex);
T.amplitude_change_percent = 100 * ...
    (T.mean_amplitude_mm - referenceAmplitude) / max(abs(referenceAmplitude), eps);
T.frequency_change_percent = 100 * ...
    (T.mean_frequency_hz - referenceFrequency) / max(abs(referenceFrequency), eps);
T.rmse_change_percent = 100 * ...
    (T.mean_plain_rmse_mV - referenceRmse) / max(abs(referenceRmse), eps);
T.reference_case = repmat(T.case_name(referenceIndex), height(T), 1);
Comparison = T;

[folder, name] = fileparts(summaryFile);
groupTag = 'all';
if nargin >= 2 && ~isempty(groupName)
    groupTag = regexprep(char(string(groupName)), '[^A-Za-z0-9_-]', '_');
end
outCsv = fullfile(folder, sprintf('%s_Comparison_%s.csv', name, groupTag));
writetable(Comparison, outCsv);
plot_comparison_local(Comparison, fullfile(folder, ...
    sprintf('%s_Comparison_%s.png', name, groupTag)));
fprintf('Reference case: %s\n', Comparison.reference_case(1));
fprintf('Saved comparison: %s\n', outCsv);
disp(Comparison(:, {'case_name','EO12_count','mean_frequency_hz', ...
    'mean_amplitude_mm','mean_plain_rmse_mV','mean_abs_delta_gap_mm', ...
    'amplitude_change_percent','frequency_change_percent', ...
    'rmse_change_percent'}));
end


function index = find_reference_local(T)
isFormal = T.window_laps == 3 & T.step_laps == 1 & ...
    T.frequency_min_hz == 300 & T.frequency_max_hz == 1000 & ...
    T.vp_top_k == 3 & abs(T.refine_half_width_hz - 2) < 1e-12 & ...
    abs(T.delta_gap_limit_mm - 0.25) < 1e-12;
index = find(isFormal, 1, 'first');
if isempty(index)
    index = 1;
end
end


function plot_comparison_local(T, outputFile)
labels = categorical(T.case_name, T.case_name, 'Ordinal', true);
fig = figure('Color', 'w', 'Visible', 'off', 'Units', 'centimeters', ...
    'Position', [2 2 23 14]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
bar(ax, labels, T.EO12_count);
ylabel(ax, 'EO12 windows'); title(ax, 'Order stability'); grid(ax, 'on');

ax = nexttile;
bar(ax, labels, T.mean_amplitude_mm);
ylabel(ax, 'Amplitude (mm)'); title(ax, 'Identified amplitude'); grid(ax, 'on');

ax = nexttile;
bar(ax, labels, T.mean_plain_rmse_mV);
ylabel(ax, 'Plain RMSE (mV)'); title(ax, 'Forward-fit residual'); grid(ax, 'on');

ax = nexttile;
bar(ax, labels, T.mean_abs_delta_gap_mm);
ylabel(ax, 'Mean |gap increment| (mm)'); title(ax, 'Gap correction'); grid(ax, 'on');

for ax = findall(fig, 'Type', 'axes').'
    ax.XTickLabelRotation = 25;
end
exportgraphics(fig, outputFile, 'Resolution', 220);
close(fig);
end
