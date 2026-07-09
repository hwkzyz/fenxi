%% Step06J_Build_FE632_TransferPrior_20241106.m
% Extract the 632 Hz FE strain-to-tip-displacement prior used by Step06K.
%
% The Workbench coordinate convention in the exported tables is:
%   X = thickness direction, Y = blade length direction, Z = width direction.
% The strain gauge position follows the measurement note:
%   X = 1.4 mm, Y = 7.05 mm, Z = 2.82 mm or 22.18 mm.

clc; clear; close all;

cfg = BTTDataConfig_20241106();

feRoot = 'D:\博士-国科\试验台数据\大台子\叶片模型';
dispFile = fullfile(feRoot, '最新位移和应变', '叶片位移(1).xls');
strainFile = fullfile(feRoot, '最新位移和应变', '叶片应变(1).xls');

outDir = fullfile(cfg.output_root, 'step06j_fe632_transfer_prior');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
figDir = fullfile(cfg.figure_root, 'step06j_fe632_transfer_prior');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

zValuesMm = [0, 1, 2.82, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.18, 24, 25].';
gaugeXmm = 1.4;
gaugeYmm = 7.05;
surfaceTolMm = 1e-6;
idwK = 8;
modeFrequencyHz = 632.3882325418;

Disp = read_fe_export_table_20241106_local(dispFile);
Strain = read_fe_export_table_20241106_local(strainFile);

[uTipMm, tipNode, tipYmm] = fe_tip_displacement_20241106_local(Disp);
uTipM = uTipMm / 1000;

Surface = Strain(abs(Strain.x_mm - gaugeXmm) <= surfaceTolMm, :);
if isempty(Surface)
    error('No FE strain surface nodes found at X = %.6g mm.', gaugeXmm);
end

nZ = numel(zValuesMm);
nearestNode = nan(nZ, 1);
nearestXmm = nan(nZ, 1);
nearestYmm = nan(nZ, 1);
nearestZmm = nan(nZ, 1);
nearestDistanceMm = nan(nZ, 1);
nearestStrain = nan(nZ, 1);
idwStrain = nan(nZ, 1);
kStrainPerMNearest = nan(nZ, 1);
kStrainPerMIdw = nan(nZ, 1);
KUmPerMicrostrainNearest = nan(nZ, 1);
KUmPerMicrostrainIdw = nan(nZ, 1);

for i = 1:nZ
    point = [gaugeXmm, gaugeYmm, zValuesMm(i)];
    info = fe_surface_idw_20241106_local(Surface, point, idwK);
    kStrainPerMNearest(i) = abs(info.nearest_strain / uTipM);
    kStrainPerMIdw(i) = abs(info.idw_strain / uTipM);

    nearestNode(i) = info.nearest_node;
    nearestXmm(i) = info.nearest_xyz(1);
    nearestYmm(i) = info.nearest_xyz(2);
    nearestZmm(i) = info.nearest_xyz(3);
    nearestDistanceMm(i) = info.nearest_distance_mm;
    nearestStrain(i) = info.nearest_strain;
    idwStrain(i) = info.idw_strain;
    KUmPerMicrostrainNearest(i) = 1 / kStrainPerMNearest(i);
    KUmPerMicrostrainIdw(i) = 1 / kStrainPerMIdw(i);
end

Rows = table(zValuesMm, nearestNode, nearestXmm, nearestYmm, nearestZmm, ...
    nearestDistanceMm, nearestStrain, idwStrain, kStrainPerMNearest, ...
    kStrainPerMIdw, KUmPerMicrostrainNearest, KUmPerMicrostrainIdw, ...
    'VariableNames', {'z_mm', 'nearest_node', 'nearest_x_mm', ...
    'nearest_y_mm', 'nearest_z_mm', 'nearest_distance_mm', ...
    'nearest_strain', 'idw_strain', 'k_strain_per_m_nearest', ...
    'k_strain_per_m_idw', 'K_um_per_microstrain_nearest', ...
    'K_um_per_microstrain_idw'});

isGaugeSide = abs(Rows.z_mm - 2.82) < 1e-9 | abs(Rows.z_mm - 22.18) < 1e-9;
KWidth = Rows.K_um_per_microstrain_idw;
KGauge = mean(Rows.K_um_per_microstrain_idw(isGaugeSide), 'omitnan');

Summary = table();
Summary.FE_source = "632.388 Hz Workbench mode-shape export";
Summary.displacement_file = string(dispFile);
Summary.strain_file = string(strainFile);
Summary.mode_frequency_hz = modeFrequencyHz;
Summary.tip_y_mm = tipYmm;
Summary.tip_node = tipNode;
Summary.u_tip_max_abs_mm = uTipMm;
Summary.gauge_x_mm = gaugeXmm;
Summary.gauge_y_mm = gaugeYmm;
Summary.gauge_z_nominal_pair_mm = "2.82 / 22.18";
Summary.FE_K_min_um_per_microstrain = min(KWidth);
Summary.FE_K_nominal_um_per_microstrain = KGauge;
Summary.FE_K_max_um_per_microstrain = max(KWidth);
Summary.FE_k_strain_per_m_at_Z2p82 = ...
    Rows.k_strain_per_m_idw(abs(Rows.z_mm - 2.82) < 1e-9);
Summary.FE_k_strain_per_m_at_Z22p18 = ...
    Rows.k_strain_per_m_idw(abs(Rows.z_mm - 22.18) < 1e-9);

summaryCsv = fullfile(outDir, 'Step06J_FE632_TransferPrior_Summary_20241106.csv');
widthCsv = fullfile(outDir, 'Step06J_FE632_TransferPrior_WidthSensitivity_20241106.csv');
matFile = fullfile(outDir, 'Step06J_FE632_TransferPrior_20241106.mat');
figFile = fullfile(figDir, 'Step06J_FE632_TransferPrior_20241106.png');

writetable(Summary, summaryCsv);
writetable(Rows, widthCsv);
save(matFile, 'Summary', 'Rows', 'Disp', 'Strain', 'dispFile', ...
    'strainFile', '-v7.3');
plot_fe_prior_20241106_local(Rows, Summary, figFile);

fprintf('\nStep06J FE 632 Hz transfer prior:\n');
disp(Summary);
fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', summaryCsv, widthCsv, ...
    matFile, figFile);


function T = read_fe_export_table_20241106_local(fileName)
opts = detectImportOptions(fileName, 'FileType', 'text', ...
    'Delimiter', '\t', 'VariableNamingRule', 'preserve');
T0 = readtable(fileName, opts);
T = table();
T.node_id = T0{:, 1};
T.x_mm = T0{:, 2};
T.y_mm = T0{:, 3};
T.z_mm = T0{:, 4};
T.value = T0{:, 5};
end


function [uTipMm, tipNode, tipYmm] = fe_tip_displacement_20241106_local(T)
tipYmm = max(T.y_mm);
tipRows = T(abs(T.y_mm - tipYmm) < 1e-9, :);
[uTipMm, ix] = max(abs(tipRows.value));
tipNode = tipRows.node_id(ix);
end


function info = fe_surface_idw_20241106_local(T, point, k)
xyz = [T.x_mm, T.y_mm, T.z_mm];
d = sqrt(sum((xyz - point).^2, 2));
[dSort, ix] = sort(d, 'ascend');
ix = ix(1:min(k, numel(ix)));
w = 1 ./ max(dSort(1:numel(ix)), 1e-9).^2;
v = T.value(ix);
idw = sum(w .* v) / sum(w);

info = struct();
info.nearest_distance_mm = dSort(1);
info.nearest_node = T.node_id(ix(1));
info.nearest_xyz = xyz(ix(1), :);
info.nearest_strain = T.value(ix(1));
info.idw_strain = idw;
end


function plot_fe_prior_20241106_local(Rows, Summary, figFile)
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 22, 12], 'Name', 'Step06J FE 632 transfer prior');
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

blue = [0.000 0.447 0.741];
red = [0.635 0.078 0.184];
gray = [0.45 0.45 0.45];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
plot(ax1, Rows.z_mm, Rows.K_um_per_microstrain_idw, '-o', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.2, ...
    'DisplayName', 'IDW');
plot(ax1, Rows.z_mm, Rows.K_um_per_microstrain_nearest, '--s', ...
    'Color', gray, 'MarkerFaceColor', 'w', 'LineWidth', 0.9, ...
    'DisplayName', 'nearest');
yline(ax1, Summary.FE_K_nominal_um_per_microstrain(1), ':', ...
    'Color', red, 'LineWidth', 1.0, 'DisplayName', 'nominal gauge pair');
xline(ax1, 2.82, ':', 'Color', red, 'HandleVisibility', 'off');
xline(ax1, 22.18, ':', 'Color', red, 'HandleVisibility', 'off');
xlabel(ax1, 'Width position Z (mm)');
ylabel(ax1, 'K (um/microstrain)');
title(ax1, 'Width-position sensitivity');
legend(ax1, 'Location', 'best');
grid(ax1, 'on');

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
plot(ax2, Rows.z_mm, abs(Rows.idw_strain), '-o', ...
    'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.2, ...
    'DisplayName', '|strain|');
yline(ax2, Summary.u_tip_max_abs_mm(1) / 1000, '--', ...
    'Color', gray, 'DisplayName', 'tip disp. (m)');
xline(ax2, 2.82, ':', 'Color', red, 'HandleVisibility', 'off');
xline(ax2, 22.18, ':', 'Color', red, 'HandleVisibility', 'off');
xlabel(ax2, 'Width position Z (mm)');
ylabel(ax2, 'FE normalized quantity');
title(ax2, 'Same mode-shape normalization');
legend(ax2, 'Location', 'best');
grid(ax2, 'on');

sgtitle(fig, sprintf(['632 Hz FE prior: K = %.3f (%.3f-%.3f) ', ...
    'um/microstrain'], Summary.FE_K_nominal_um_per_microstrain(1), ...
    Summary.FE_K_min_um_per_microstrain(1), ...
    Summary.FE_K_max_um_per_microstrain(1)));
set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 9);
exportgraphics(fig, figFile, 'Resolution', 300);
close(fig);
end
