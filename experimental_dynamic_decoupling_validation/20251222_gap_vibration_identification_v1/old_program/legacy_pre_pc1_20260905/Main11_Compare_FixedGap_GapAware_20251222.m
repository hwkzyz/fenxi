%% Main11: Compare independent fixed-gap and gap-aware identification
clear; clc; close all;
cfg = Setup_Paths_20251222();
C = CaseConfig();
sensorTag = ['S', sprintf('%d', C.sensorIds)];
regionTag = lower(C.flowConfig.identification.resonanceRegionShortTag);

fixedRegionRoot = fullfile(cfg.paths.foundationResults, ...
    ['step05_noeta_vptopk_direct_template_', regionTag, '*']);
fixedFile = newest_result_local(fixedRegionRoot, ...
    sprintf('Result_Step05_NoEtaVPTopKDirectTemplate_B%d_%s_20251222.mat', ...
    C.bladeId, sensorTag));
gapFile = newest_result_local(cfg.paths.gapResults, ...
    sprintf('Main_GapAware_VPTopK_FullWave_20251222_B%d_%s_%s_*.mat', ...
    C.bladeId, sensorTag, regionTag));
assert(~isempty(fixedFile), 'Run Main09 first; no fixed-gap result was found.');
assert(~isempty(gapFile), 'Run Main10 first; no gap-aware result was found.');

F = load(fixedFile, 'Result');
G = load(gapFile, 'Result');
Tf = fixed_trend_local(F.Result);
Tg = gap_trend_local(G.Result);
n = min(height(Tf), height(Tg));
Tf = Tf(1:n, :);
Tg = Tg(1:n, :);
Comparison = table(Tf.Window, Tf.TimeSec, ...
    Tf.EO, Tf.FrequencyHz, Tf.AmplitudeMm, Tf.PlainRmseMv, ...
    Tg.EO, Tg.FrequencyHz, Tg.AmplitudeMm, Tg.PlainRmseMv, Tg.MeanDeltaGapMm, ...
    'VariableNames', {'Window','TimeSec','FixedEO','FixedFrequencyHz', ...
    'FixedAmplitudeMm','FixedPlainRmseMv','GapEO','GapFrequencyHz', ...
    'GapAmplitudeMm','GapPlainRmseMv','GapMeanDeltaGapMm'});

csvFile = fullfile(cfg.paths.comparison, sprintf( ...
    'Comparison_20251222_B%d_%s_%s.csv', C.bladeId, sensorTag, regionTag));
matFile = fullfile(cfg.paths.comparison, sprintf( ...
    'Comparison_20251222_B%d_%s_%s.mat', C.bladeId, sensorTag, regionTag));
writetable(Comparison, csvFile);
save(matFile, 'Comparison', 'fixedFile', 'gapFile');
fprintf('Compared %d independent windows.\nFixed: %s\nGap:   %s\nCSV:   %s\n', ...
    n, fixedFile, gapFile, csvFile);

function file = newest_result_local(rootDir, namePattern)
files = dir(fullfile(rootDir, '**', namePattern));
if isempty(files), file = ''; return; end
[~, idx] = max([files.datenum]);
file = fullfile(files(idx).folder, files(idx).name);
end

function T = fixed_trend_local(R)
W = R.WindowResult(:); n = numel(W);
window = (1:n).'; time = nan(n,1); eo = nan(n,1); f = nan(n,1);
a = nan(n,1); prmse = nan(n,1);
for i = 1:n
    q = W(i).Result;
    if isfield(W(i),'time_window_s'), time(i)=mean(W(i).time_window_s); end
    eo(i)=q.EO_id; f(i)=q.fn_id; a(i)=q.A_id;
    if isfield(q,'plain_voltage_rmse'), prmse(i)=1000*q.plain_voltage_rmse; end
end
T = table(window,time,eo,f,a,prmse,'VariableNames', ...
    {'Window','TimeSec','EO','FrequencyHz','AmplitudeMm','PlainRmseMv'});
end

function T = gap_trend_local(R)
W = R.WindowResult(:); n = numel(W);
window = (1:n).'; time = nan(n,1); eo = nan(n,1); f = nan(n,1);
a = nan(n,1); prmse = nan(n,1); dg = nan(n,1);
for i = 1:n
    q = W(i).modelFits.gap_only;
    if isfield(W(i),'windowId'), window(i)=W(i).windowId; end
    if isfield(W(i),'timeWindow'), time(i)=mean(W(i).timeWindow); end
    eo(i)=q.EO; f(i)=q.freqHz; a(i)=q.amplitudeMm; prmse(i)=q.plainRmseMv;
    if isfield(q,'deltaGapMm'), dg(i)=mean(q.deltaGapMm,'omitnan'); end
end
T = table(window,time,eo,f,a,prmse,dg,'VariableNames', ...
    {'Window','TimeSec','EO','FrequencyHz','AmplitudeMm','PlainRmseMv','MeanDeltaGapMm'});
end
