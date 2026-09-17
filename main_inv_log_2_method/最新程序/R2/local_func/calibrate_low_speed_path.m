function [pathCal, fit] = calibrate_low_speed_path(templateLib, low, opts)
%CALIBRATE_LOW_SPEED_PATH  Fit a bounded low-speed clearance path.
if nargin < 3 || isempty(opts), opts = struct(); end
xAll = low.xGrid(:);
yAll = low.templateLow(:);
if isfield(opts, 'xDomain') && numel(opts.xDomain) == 2
    keep = xAll >= min(opts.xDomain) & xAll <= max(opts.xDomain);
else
    keep = true(size(xAll));
end
x = xAll(keep);
y = yAll(keep);
if isfield(low, 'stdLow') && isfield(low, 'countLow')
    sigma = low.stdLow(keep) ./ sqrt(max(low.countLow(keep), 1));
else
    sigma = ones(size(y));
end
sigmaFloor = get_opt(opts, 'sigmaFloor', max(std(y) * 1e-3, eps));
weights = 1 ./ max(sigma, sigmaFloor);

gTrain = templateLib.gapTrain(:);
gMin = min(gTrain);
gMax = max(gTrain);
theta0 = [get_opt(opts, 'g0Init', mean([gMin, gMax])), ...
    get_opt(opts, 'muInit', 0), get_opt(opts, 'tauInit', 0), ...
    get_opt(opts, 'zetaInit', 1)];
lb = [get_opt(opts, 'g0Lower', gMin), get_opt(opts, 'muLower', -0.15), ...
    get_opt(opts, 'tauLower', min(x)), get_opt(opts, 'zetaLower', 0.8)];
ub = [get_opt(opts, 'g0Upper', gMax), get_opt(opts, 'muUpper', 0.15), ...
    get_opt(opts, 'tauUpper', max(x)), get_opt(opts, 'zetaUpper', 1.2)];

% Some validation cases have a known geometric path.  In that case only the
% requested path parameters are optimized; the remaining ones stay fixed.
% This prevents a no-tilt calibration from explaining response-model mismatch
% by inventing a local gradient or a spatial re-registration.
fitMask = logical(get_opt(opts, 'fitMask', true(1, 4)));
if numel(fitMask) ~= 4
    error('opts.fitMask must contain four logical values for [g0 mu tau zeta].');
end
theta = theta0;
if any(fitMask)
    free0 = theta0(fitMask);
    freeLb = lb(fitMask);
    freeUb = ub(fitMask);
    resFun = @(freeTheta) weighted_path_residual_masked(freeTheta, theta0, ...
        fitMask, templateLib, x, y, weights, get_opt(opts, 'fitGainOffset', true));
    [freeTheta, solveInfo] = solve_lsq_bounded(resFun, free0, freeLb, freeUb, ...
        get_opt(opts, 'maxIter', 500));
    theta(fitMask) = freeTheta;
else
    solveInfo = struct('exitflag', 0, 'iterations', 0, ...
        'message', 'All path parameters fixed by opts.fitMask.');
end
rPath = eval_gap_path_response(templateLib, theta(1), theta(2), ...
    theta(3), theta(4), x);
if get_opt(opts, 'fitGainOffset', true)
    [b, kappa] = fit_gain_offset(rPath, y, weights);
else
    b = get_opt(opts, 'baseline', 0);
    kappa = get_opt(opts, 'kappa', 1);
end
residual = weights .* (y - (b + kappa .* rPath));

pathCal = struct('g0', theta(1), 'mu', theta(2), 'tau', theta(3), ...
    'zeta', theta(4), 'kappa', kappa, 'b', b, ...
    'calibration_source', "low_speed_template");
fit = struct('x', x, 'y', y, 'path_response', rPath, ...
    'fitted', b + kappa .* rPath, 'residual', residual, ...
    'rmse_V', sqrt(mean((y - (b + kappa .* rPath)).^2)), ...
    'weighted_rmse', sqrt(mean(residual.^2)), ...
    'solve_info', solveInfo, 'theta', theta);
end

function r = weighted_path_residual_masked(freeTheta, thetaFixed, fitMask, ...
    templateLib, x, y, weights, fitGainOffset)
theta = thetaFixed;
theta(fitMask) = freeTheta;
r = weighted_path_residual(theta, templateLib, x, y, weights, fitGainOffset);
end

function r = weighted_path_residual(theta, templateLib, x, y, weights, fitGainOffset)
rPath = eval_gap_path_response(templateLib, theta(1), theta(2), ...
    theta(3), theta(4), x);
if fitGainOffset
    [b, kappa] = fit_gain_offset(rPath, y, weights);
else
    b = 0;
    kappa = 1;
end
r = weights .* (y - (b + kappa .* rPath));
end

function [b, kappa] = fit_gain_offset(rPath, y, weights)
X = [ones(size(rPath)), rPath];
W = weights(:);
normal = (X .* W)' * (X .* W);
rhs = (X .* W)' * (y(:) .* W);
coef = normal \ rhs;
b = coef(1);
kappa = max(coef(2), eps);
end

function value = get_opt(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end
