function pathModel = build_path_template_model(templateLib, templateLow, pathCal)
%BUILD_PATH_TEMPLATE_MODEL  Package a fixed low-speed clearance path.
%
% pathCal must contain g0, mu, tau, zeta, and kappa.  The path parameters
% are calibration quantities and are not optimized by the high-speed solver.

required = {'g0', 'mu', 'tau', 'zeta', 'kappa'};
for i = 1:numel(required)
    if ~isfield(pathCal, required{i}) || ~isscalar(pathCal.(required{i})) || ...
            ~isfinite(pathCal.(required{i}))
        error('pathCal.%s must be a finite scalar.', required{i});
    end
end
if pathCal.zeta <= 0
    error('pathCal.zeta must be positive.');
end

xTemplate = templateLib.xGrid(:);
if isvector(templateLow), templateLow = templateLow(:); end
if numel(xTemplate) < 2 || numel(unique(xTemplate)) < 2
    error('templateLib.xGrid must contain at least two distinct spatial samples.');
end
if ~isfield(templateLib, 'gapTrain') || numel(unique(templateLib.gapTrain(:))) < 2
    error('templateLib.gapTrain must contain at least two distinct gap samples.');
end
if size(templateLow, 1) ~= numel(xTemplate)
    error('templateLow rows must match templateLib.xGrid.');
end
if any(~isfinite(templateLow))
    error('templateLow contains non-finite values.');
end

pathModel = struct();
pathModel.templateLib = templateLib;
pathModel.baseLib = templateLib;
pathModel.domain = templateLib.domain;
pathModel.xGrid = templateLib.xGrid;
pathModel.gapTrain = templateLib.gapTrain;
pathModel.xTemplate = xTemplate;
pathModel.templateLow = templateLow;
% An adaptive spline estimator can supply a stable low-speed derivative.
% Keeping it alongside the sampled forward template lets EO screening use
% the continuous derivative without changing the final nonlinear forward
% waveform model.
if isfield(pathCal,'lowDerivativeBySensor') && ...
        isequal(size(pathCal.lowDerivativeBySensor),size(templateLow))
    pathModel.lowDerivativeBySensor=pathCal.lowDerivativeBySensor;
end
if size(templateLow,2)>1
    if isfield(pathCal,'sensorIds') && numel(pathCal.sensorIds)==size(templateLow,2)
        pathModel.sensorIds=pathCal.sensorIds(:).';
    else
        pathModel.sensorIds=1:size(templateLow,2);
    end
else
    pathModel.sensorIds=1;
end
pathModel.pathCal = pathCal;
pathModel.model_type = "fixed_low_speed_path_increment";
pathModel.path_fixed = true;
% Keep the public evaluator contract explicit.  eval_gap_template checks
% this flag to dispatch to the low-speed-template increment model.
pathModel.fixedPathIncrement = true;
% Cache the absolute response surface used by the high-speed evaluator. The
% direct pointwise fallback is correct but too slow for VP scans.
gLo = max(min(templateLib.gapTrain), pathCal.g0 - 0.75);
gHi = min(max(templateLib.gapTrain), pathCal.g0 + 0.75);
gCache = linspace(gLo, gHi, 161).';
if numel(unique(gCache)) < 2
    error('The fixed-path cache requires a non-degenerate gap interval.');
end
[sCache, ~] = eval_gap_grid(templateLib, gCache);
pathModel.pathCache = struct('gGrid', gCache, ...
    'xGrid', templateLib.xGrid(:), 'SGrid', sCache);
pathModel.pathCache.interpolant = griddedInterpolant( ...
    {gCache, templateLib.xGrid(:)}, sCache, 'linear', 'none');
end
