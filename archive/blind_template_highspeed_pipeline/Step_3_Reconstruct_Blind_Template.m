%% Step 3: reconstruct the shared high-speed template from local segments
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'blind_template_highspeed_results');
load(fullfile(outDir, 'stage0_config.mat'), 'cfg');
load(fullfile(outDir, 'stage2_segments.mat'), 'segments', 'xRef');

if isempty(segments)
    error('No valid local segments found. Step 2 must produce at least one segment.');
end

[T0, seedInfo] = initialize_template_from_segments(segments, cfg.blind);
[blindResult, fitList] = reconstruct_blind_template(segments, xRef, T0, cfg.blind);
blindResult.seedInfo = seedInfo;
blindResult.fitList = fitList;

save(fullfile(outDir, 'stage3_blind_template.mat'), ...
    'blindResult', 'fitList', '-v7.3');

figure('Name', 'Blind-Template Step 3 - Reconstruction', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 22, 11]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on;
plot(xRef, T0, 'c--', 'LineWidth', 1.0, 'DisplayName', 'initial template');
plot(xRef, blindResult.T, 'b-', 'LineWidth', 1.4, 'DisplayName', 'reconstructed template');
xlabel('Local coordinate x (mm)');
ylabel('Voltage');
title('Shared template reconstruction');
legend('Location', 'best', 'Box', 'off');

nexttile;
plot(blindResult.costHist, 'o-', 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('Mean registration RMSE');
title('Iteration cost history');

nexttile;
histogram(blindResult.d, 20);
xlabel('Recovered shift d_q (mm)');
ylabel('Count');
title('Recovered local shifts');

nexttile;
scatter([segments.trueShift]', blindResult.d, 18, [segments.sensor]', 'filled');
xlabel('Truth shift u(t_q) (mm)');
ylabel('Recovered shift d_q (mm)');
title('Recovered shift vs truth');
grid on;

fprintf('[Step 3] Blind-template reconstruction finished in %d iterations.\n', blindResult.numIter);
fprintf('[Step 3] Final mean registration RMSE = %.6g\n', blindResult.costHist(end));

function [T0, seedInfo] = initialize_template_from_segments(segments, blindCfg)
weights = [segments.weight]';
[~, order] = sort(weights, 'descend');
seedIds = order(1:min(blindCfg.seedCount, numel(order)));
seedTemplate = segments(seedIds(1)).yInterp;
seedMask = segments(seedIds(1)).mask;
seedTemplate = local_fill_missing(seedTemplate, seedMask);

fitSeed = repmat(struct('d', 0, 'a', 1, 'b', 0, 'rmse', inf), numel(seedIds), 1);
fitSeed(1).d = 0;
fitSeed(1).a = 1;
fitSeed(1).b = 0;
fitSeed(1).rmse = 0;
for i = 2:numel(seedIds)
    fitSeed(i) = register_segment_to_template(segments(seedIds(i)), seedTemplate, blindCfg);
end
T0 = update_template_from_segments(segments(seedIds), fitSeed, blindCfg);
seedInfo = struct('seedIds', seedIds(:), 'seedPrimary', seedIds(1));
end

function [result, fitList] = reconstruct_blind_template(segments, xRef, T0, blindCfg)
T = T0(:);
costHist = zeros(blindCfg.maxIter, 1);
templateChangeHist = zeros(blindCfg.maxIter, 1);
fitList = repmat(struct('d', 0, 'a', 1, 'b', 0, 'rmse', inf), numel(segments), 1);

for iter = 1:blindCfg.maxIter
    for q = 1:numel(segments)
        fitList(q) = register_segment_to_template(segments(q), T, blindCfg);
    end

    dVec = [fitList.d]';
    dMean = mean(dVec);
    for q = 1:numel(fitList)
        fitList(q).d = fitList(q).d - dMean;
    end

    Tnew = update_template_from_segments(segments, fitList, blindCfg);
    templateChangeHist(iter) = norm(Tnew - T) / max(norm(T), eps);
    costHist(iter) = mean([fitList.rmse]);
    T = Tnew;
    if templateChangeHist(iter) < blindCfg.tolTemplate
        costHist = costHist(1:iter);
        templateChangeHist = templateChangeHist(1:iter);
        break;
    end
end

result = struct();
result.T = T(:);
result.xRef = xRef(:);
result.d = [fitList.d]';
result.a = [fitList.a]';
result.b = [fitList.b]';
result.rmse = [fitList.rmse]';
result.costHist = costHist(:);
result.templateChangeHist = templateChangeHist(:);
result.numIter = numel(costHist);
end

function fit = register_segment_to_template(seg, T, blindCfg)
y = seg.yInterp(:);
mask = seg.mask(:) & isfinite(T(:)) & isfinite(y);
T = T(:);

best = struct('d', 0, 'a', 1, 'b', 0, 'rmse', inf);
dGridCoarse = blindCfg.dGridCoarse(:)';
for d = dGridCoarse
    trial = solve_shift_gain_bias(seg.xRef, y, mask, T, d, blindCfg.useGainBias);
    if trial.rmse < best.rmse
        best = trial;
    end
end

dFineLo = best.d - blindCfg.dFineHalfWidth;
dFineHi = best.d + blindCfg.dFineHalfWidth;
dGridFine = dFineLo:blindCfg.dFineStep:dFineHi;
for d = dGridFine
    trial = solve_shift_gain_bias(seg.xRef, y, mask, T, d, blindCfg.useGainBias);
    if trial.rmse < best.rmse
        best = trial;
    end
end
fit = best;
end

function trial = solve_shift_gain_bias(xRef, y, mask, T, d, useGainBias)
Tshift = interp1(xRef - d, T, xRef, 'pchip', NaN);
valid = mask & isfinite(Tshift);
if nnz(valid) < 12
    trial = struct('d', d, 'a', 1, 'b', 0, 'rmse', inf);
    return;
end
if useGainBias
    A = [Tshift(valid), ones(nnz(valid), 1)];
    theta = A \ y(valid);
    a = theta(1);
    b = theta(2);
else
    a = 1;
    b = 0;
end
res = y(valid) - (a * Tshift(valid) + b);
trial = struct('d', d, 'a', a, 'b', b, 'rmse', sqrt(mean(res.^2)));
end

function Tnew = update_template_from_segments(segments, fitList, blindCfg)
xRef = segments(1).xRef(:);
num = zeros(size(xRef));
den = zeros(size(xRef));
for q = 1:numel(segments)
    y = segments(q).yInterp(:);
    mask = segments(q).mask(:);
    d = fitList(q).d;
    a = fitList(q).a;
    b = fitList(q).b;
    if abs(a) < 1e-8
        continue;
    end
    xValid = xRef(mask);
    yValid = y(mask);
    yWarp = interp1(xValid + d, yValid, xRef, 'pchip', NaN);
    valid = isfinite(yWarp);
    if nnz(valid) < 12
        continue;
    end
    yNorm = (yWarp(valid) - b) / a;
    w = segments(q).weight;
    num(valid) = num(valid) + w * yNorm;
    den(valid) = den(valid) + w;
end
Tnew = num ./ max(den, eps);
missing = ~isfinite(Tnew) | den <= 0;
Tnew = local_fill_missing(Tnew, ~missing);
Tnew = smoothdata(Tnew, 'sgolay', blindCfg.smoothWindow);
end

function y = local_fill_missing(y, mask)
x = (1:numel(y))';
valid = mask(:) & isfinite(y(:));
if nnz(valid) < 2
    y(~valid) = 0;
    return;
end
y(~valid) = interp1(x(valid), y(valid), x(~valid), 'linear', 'extrap');
end
