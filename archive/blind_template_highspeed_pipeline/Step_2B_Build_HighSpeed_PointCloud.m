%% Step 2B: keep all high-speed passing samples as a point cloud
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage1_simulated_data.mat'), 'Data_High', 'truth');

[cloud, activeThreshold] = build_highspeed_point_cloud(Data_High, truth, cfg);
save(fullfile(outDir, 'stage2B_point_cloud.mat'), 'cloud', 'activeThreshold', '-v7.3');

if cfg.makeFigures
    figure('Name', 'Blind-Template Step 2B - High-Speed Point Cloud', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, 22, 10]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    scatter(cloud.x, cloud.V, 6, cloud.t * 1000, 'filled');
    xlabel('Nominal coordinate x (mm)');
    ylabel('Voltage');
    title('Mapped high-speed point cloud');
    colorbar;

    nexttile;
    histogram(cloud.x, 80);
    xlabel('Nominal coordinate x (mm)');
    ylabel('Count');
    title(sprintf('All active high-speed samples, N = %d', numel(cloud.V)));
end

fprintf('[Step 2B] Point-cloud samples: %d\n', numel(cloud.V));
fprintf('[Step 2B] Passing events represented: %d\n', numel(unique([cloud.rev, cloud.sensor], 'rows')));

function [cloud, activeThreshold] = build_highspeed_point_cloud(Data_High, truth, cfg)
T_opr = Data_High.T_opr_truth(:);
t = Data_High.t(:);
V = Data_High.V_cap(:);
xHalf = 0.50 * diff(truth.domain);
base = local_percentile(V, 5);
span = local_percentile(V, 99) - base;
activeThreshold = base + cfg.segment.activeLevel * span;

tAll = [];
xAll = [];
vAll = [];
revAll = [];
sensorAll = [];
eventAll = [];
eventId = 0;

for m = 1:numel(T_opr)-1
    Omega = 2*pi / (T_opr(m+1) - T_opr(m));
    for s = 1:numel(cfg.alpha_k)
        eventId = eventId + 1;
        Te = T_opr(m) + cfg.alpha_k(s) / Omega;
        idx = find(abs(t - Te) <= xHalf / (Omega * cfg.R_tip));
        if isempty(idx), continue; end
        xLocal = Omega * cfg.R_tip * (t(idx) - Te);
        keep = xLocal >= truth.domain(1) & xLocal <= truth.domain(2) & V(idx) > activeThreshold;
        idx = idx(keep);
        xLocal = xLocal(keep);
        if numel(idx) < cfg.segment.minPointCount
            continue;
        end
        tAll = [tAll; t(idx)]; %#ok<AGROW>
        xAll = [xAll; xLocal(:)]; %#ok<AGROW>
        vAll = [vAll; V(idx)]; %#ok<AGROW>
        revAll = [revAll; m * ones(numel(idx), 1)]; %#ok<AGROW>
        sensorAll = [sensorAll; s * ones(numel(idx), 1)]; %#ok<AGROW>
        eventAll = [eventAll; eventId * ones(numel(idx), 1)]; %#ok<AGROW>
    end
end

xGrid = linspace(truth.domain(1), truth.domain(2), cfg.profile.xGridN)';
cloud = struct();
cloud.t = tAll(:);
cloud.x = xAll(:);
cloud.V = vAll(:);
cloud.rev = revAll(:);
cloud.sensor = sensorAll(:);
cloud.event = eventAll(:);
cloud.domain = truth.domain;
cloud.xGrid = xGrid;
end

function p = local_percentile(x, pct)
x = sort(x(isfinite(x)));
if isempty(x), p = NaN; return; end
q = 1 + (numel(x) - 1) * pct / 100;
lo = floor(q);
hi = ceil(q);
if lo == hi
    p = x(lo);
else
    p = x(lo) + (q - lo) * (x(hi) - x(lo));
end
end
