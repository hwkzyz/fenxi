%% Plot 05: fixed-trust static and projected gap-cost curves.
close all;
clc;

mainDir = fileparts(fileparts(mfilename('fullpath')));
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);

if ~exist('gTrue', 'var') || isempty(gTrue)
    gTrue = 0.2;
end
D = compute_inv_log_2_fixed_demo(ctx, gTrue);

gapGrid = linspace(max(0.05, D.staticState.gHat - 0.12), ...
    min(max(D.templateInv.gapTrain) + D.cfgCase.rawGapSearchMargin, ...
    D.staticState.gHat + 0.12), 91);
[JStatic, JProj] = compute_cost_curves_local(D.highMap, D.templateInv, ...
    D.staticState, D.gapState.freq_list(:), gapGrid);

fig = figure('Name', 'fixed trust gap cost curves', 'Color', 'w');
ax = axes(fig);
plot(ax, gapGrid, JStatic, '--', 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 1.2, 'DisplayName', 'static cost J_{static}(g)');
hold(ax, 'on');
plot(ax, gapGrid, JProj, '-', 'Color', [0.00 0.45 0.74], ...
    'LineWidth', 1.5, 'DisplayName', 'projected cost J_{proj}(g)');
xline(ax, gTrue, 'k--', 'LineWidth', 0.9, 'DisplayName', 'true g');
xline(ax, D.staticState.gHat, '-', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 1.1, 'DisplayName', 'static estimate');
xline(ax, D.gapState.gHat, '-', 'Color', [0.47 0.67 0.19], ...
    'LineWidth', 1.1, 'DisplayName', 'projected estimate');
xlabel(ax, 'Candidate gap g (mm)', 'FontSize', 9);
ylabel(ax, 'Residual cost', 'FontSize', 9);
title(ax, 'Gap calculation with fixed trusted window', 'FontWeight', 'normal', 'FontSize', 9);
legend(ax, 'Location', 'best');
apply_inv_log_2_figure_style(ax, 170, 86);

function [JStatic, JProj] = compute_cost_curves_local(highMap, templateLib, staticState, freqList, gapGrid)
x = highMap.x_v(:);
t = highMap.t_v(:);
V = highMap.V_a(:);
dx0 = staticState.dx0;
JStatic = zeros(numel(gapGrid), 1);
JProj = zeros(numel(gapGrid), 1);
for ig = 1:numel(gapGrid)
    g = gapGrid(ig);
    F = eval_gap_template(templateLib, g, x - dx0);
    valid = isfinite(F) & isfinite(V);
    r = V(valid) - F(valid);
    JStatic(ig) = mean(r .^ 2);
    Fx = eval_gap_derivative(templateLib, g, x(valid) - dx0);
    B = build_basis_local(Fx, t(valid), freqList);
    coef = B \ r;
    rProj = r - B * coef;
    JProj(ig) = mean(rProj .^ 2);
end
end

function B = build_basis_local(Fx, t, freqList)
cols = {Fx};
for ifr = 1:numel(freqList)
    w = 2 * pi * freqList(ifr);
    cols{end + 1} = Fx .* sin(w * t); %#ok<AGROW>
    cols{end + 1} = Fx .* cos(w * t); %#ok<AGROW>
end
B = [cols{:}];
end
