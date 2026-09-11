%% Main11: compare the formal fixed-gap and gap-aware results.
clear; clc; close all;

cfg = Setup_Paths_20250527();
foundationFile = strtrim(getenv('FORMAL_COMPARE_FOUNDATION_FILE'));
if isempty(foundationFile)
    foundationFile = cfg.files.formalFoundationResult;
end
gapFile = strtrim(getenv('FORMAL_COMPARE_GAP_FILE'));
if isempty(gapFile)
    gapFile = cfg.files.formalGapResult;
end
comparisonTag = regexprep(strtrim(getenv('FORMAL_COMPARE_TAG')), '[^\w\-]', '_');
assert(isfile(foundationFile), 'Missing formal Foundation result: %s', foundationFile);
assert(isfile(gapFile), 'Missing formal GapAware result: %s', gapFile);

F = load(foundationFile, 'Result');
G = load(gapFile, 'Result');
assert(isfield(F, 'Result') && isfield(F.Result, 'Trend'), ...
    'Foundation result does not contain Result.Trend.');
assert(isfield(G, 'Result') && isfield(G.Result, 'Trend'), ...
    'GapAware result does not contain Result.Trend.');

TF = F.Result.Trend;
TG = G.Result.Trend;
requiredFoundation = {'window_id','EO_id','fn_id','A_id','phi_id_wrapped', ...
    'dx_c_id','plain_voltage_rmse'};
% Historical Step07J exports used the name gap_waveform_rmse_mV;
% accept that frozen equivalent and normalize it below.
if ~ismember('gap_plain_rmse_mV', TG.Properties.VariableNames) && ...
        ismember('gap_waveform_rmse_mV', TG.Properties.VariableNames)
    TG.gap_plain_rmse_mV = TG.gap_waveform_rmse_mV;
end
requiredGap = {'window_id','gap_EO','gap_frequency_hz','gap_amplitude_mm', ...
    'gap_dx_mm','gap_plain_rmse_mV','gap_mean_delta_gap_mm'};
assert(all(ismember(requiredFoundation, TF.Properties.VariableNames)), ...
    'Foundation trend lacks required formal columns.');
assert(all(ismember(requiredGap, TG.Properties.VariableNames)), ...
    'GapAware trend lacks required formal columns.');

foundation = table(TF.window_id, TF.EO_id, TF.fn_id, TF.A_id, ...
    TF.phi_id_wrapped, TF.dx_c_id, 1000 .* TF.plain_voltage_rmse, ...
    'VariableNames', {'window_id','foundation_EO','foundation_frequency_hz', ...
    'foundation_amplitude_mm','foundation_phase_rad','foundation_dx_mm', ...
    'foundation_plain_rmse_mV'});
gap = table(TG.window_id, TG.gap_EO, TG.gap_frequency_hz, ...
    TG.gap_amplitude_mm, TG.gap_dx_mm, TG.gap_plain_rmse_mV, ...
    TG.gap_mean_delta_gap_mm, ...
    'VariableNames', {'window_id','gap_EO','gap_frequency_hz', ...
    'gap_amplitude_mm','gap_dx_mm','gap_plain_rmse_mV', ...
    'gap_mean_delta_gap_mm'});
if ismember('gap_frequency_offset_hz', TG.Properties.VariableNames)
    gap.gap_frequency_offset_hz = TG.gap_frequency_offset_hz;
end

WindowComparison = innerjoin(foundation, gap, 'Keys', 'window_id');
assert(height(WindowComparison) == cfg.window.count, ...
    'Expected %d common windows, found %d.', cfg.window.count, height(WindowComparison));
WindowComparison.delta_EO = WindowComparison.gap_EO - WindowComparison.foundation_EO;
WindowComparison.delta_frequency_hz = WindowComparison.gap_frequency_hz - ...
    WindowComparison.foundation_frequency_hz;
WindowComparison.delta_amplitude_mm = WindowComparison.gap_amplitude_mm - ...
    WindowComparison.foundation_amplitude_mm;
WindowComparison.delta_plain_rmse_mV = WindowComparison.gap_plain_rmse_mV - ...
    WindowComparison.foundation_plain_rmse_mV;

Summary = table(mean(WindowComparison.foundation_amplitude_mm, 'omitnan'), ...
    mean(WindowComparison.gap_amplitude_mm, 'omitnan'), ...
    mean(WindowComparison.delta_amplitude_mm, 'omitnan'), ...
    mean(WindowComparison.foundation_plain_rmse_mV, 'omitnan'), ...
    mean(WindowComparison.gap_plain_rmse_mV, 'omitnan'), ...
    nnz(WindowComparison.foundation_EO == WindowComparison.gap_EO), ...
    height(WindowComparison), ...
    'VariableNames', {'foundation_mean_amplitude_mm','gap_mean_amplitude_mm', ...
    'mean_amplitude_change_mm','foundation_mean_plain_rmse_mV', ...
    'gap_mean_plain_rmse_mV','same_EO_window_count','window_count'});

if ~isfolder(cfg.paths.comparison)
    mkdir(cfg.paths.comparison);
end
if isempty(comparisonTag)
    comparisonStem = 'Formal_Comparison_20250527_B1_S136';
else
    comparisonStem = ['Formal_Comparison_20250527_B1_S136_', comparisonTag];
end
matFile = fullfile(cfg.paths.comparison, [comparisonStem, '.mat']);
csvFile = fullfile(cfg.paths.comparison, [comparisonStem, '_Windows.csv']);
summaryFile = fullfile(cfg.paths.comparison, [comparisonStem, '_Summary.csv']);
save(matFile, 'WindowComparison', 'Summary', 'foundationFile', 'gapFile');
writetable(WindowComparison, csvFile);
writetable(Summary, summaryFile);

if cfg.run.showPlots || cfg.run.saveFigures
    fig = figure('Color', 'w', 'Name', 'Formal fixed-gap and gap-aware comparison');
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    plot(WindowComparison.window_id, WindowComparison.foundation_amplitude_mm, '-o'); hold on;
    plot(WindowComparison.window_id, WindowComparison.gap_amplitude_mm, '-s');
    ylabel('Amplitude (mm)'); legend('Foundation', 'GapAware', 'Location', 'best'); grid on;
    nexttile;
    plot(WindowComparison.window_id, WindowComparison.foundation_frequency_hz, '-o'); hold on;
    plot(WindowComparison.window_id, WindowComparison.gap_frequency_hz, '-s');
    ylabel('Frequency (Hz)'); grid on;
    nexttile;
    plot(WindowComparison.window_id, WindowComparison.foundation_plain_rmse_mV, '-o'); hold on;
    plot(WindowComparison.window_id, WindowComparison.gap_plain_rmse_mV, '-s');
    xlabel('Window'); ylabel('Plain RMSE (mV)'); grid on;
    if cfg.run.saveFigures
        exportgraphics(fig, fullfile(cfg.paths.comparison, ...
            [comparisonStem, '.png']), 'Resolution', 300);
    end
end

disp(Summary);
fprintf('Saved formal comparison:\n  %s\n  %s\n', matFile, csvFile);
