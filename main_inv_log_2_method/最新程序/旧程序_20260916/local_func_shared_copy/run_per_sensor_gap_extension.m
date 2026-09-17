function R = run_per_sensor_gap_extension(highMap, templateLib, cfg, lowState)
%RUN_PER_SENSOR_GAP_EXTENSION Experimental per-sensor gap extension.
%
% This is deliberately not the default solver. It keeps the Route-30 shared
% VP result as a frequency seed, then fits one absolute gap per sensor in the
% complete voltage model. The extension currently targets the absolute
% template mode; low-speed path-increment sensor models are a separate stage.

if isfield(templateLib, 'fixedPathIncrement') && templateLib.fixedPathIncrement
    error('invlog2:PerSensorPathUnsupported', ...
        'Per-sensor extension currently requires absolute templates.');
end
if nargin < 4 || isempty(lowState) || ~isfield(lowState,'gHat')
    error('invlog2:MissingLowSpeedState', ...
        'Per-sensor extension requires an explicit low-speed state.');
end

shared = run_inv_log_2_main_method(highMap, templateLib, cfg, lowState);
if shared.model_order == 1
    error('invlog2:PerSensorDualRequired', ...
        'Per-sensor extension currently requires a dual-frequency record.');
end
[fitMap, queryInfo] = apply_query_safe_window_mask(highMap, templateLib, cfg);
t = fitMap.t_v(:); x = fitMap.x_v(:); V = fitMap.V_a(:);
sensorValues = unique(fitMap.S_v(:), 'stable');
S = numel(sensorValues);
if S < 2
    error('invlog2:PerSensorInsufficientSensors', ...
        'At least two sensors are required for the extension.');
end
[~, sensorId] = ismember(fitMap.S_v(:), sensorValues);

fSeed = shared.f_id(:).';
pSeed = shared.fit.p(:).';
if isfield(cfg,'perSensorFrequencySeed') && numel(cfg.perSensorFrequencySeed)==2
    fSeed = cfg.perSensorFrequencySeed(:).';
    pSeed([3,6]) = fSeed;
end
gSeed = repmat(shared.g_used, 1, S);
z0 = [gSeed, shared.dx_used, pSeed];
gapHalf = get_field(cfg,'perSensorGapHalfWidthMm',0.12);
dxHalf = get_field(cfg,'perSensorDxHalfWidthMm',0.20);
ampUpper = get_field(cfg,'nonlinearAmpUpperBound',0.8);
freqHalf = get_field(cfg,'perSensorFrequencyHalfWidthHz',80);
freqMin = min([cfg.f1Grid(:);cfg.f2Grid(:)]);
freqMax = max([cfg.f1Grid(:);cfg.f2Grid(:)]);
lb = [gSeed-gapHalf, shared.dx_used-dxHalf, ...
    max(0.001,pSeed(1)-ampUpper), -pi, max(freqMin,fSeed(1)-freqHalf), ...
    max(0.001,pSeed(4)-ampUpper), -pi, max(freqMin,fSeed(2)-freqHalf)];
ub = [gSeed+gapHalf, shared.dx_used+dxHalf, ...
    min(ampUpper,pSeed(1)+ampUpper), pi, min(freqMax,fSeed(1)+freqHalf), ...
    min(ampUpper,pSeed(4)+ampUpper), pi, min(freqMax,fSeed(2)+freqHalf)];
lb(1:S) = max(lb(1:S), min(templateLib.gapTrain)-cfg.rawGapSearchMargin);
ub(1:S) = min(ub(1:S), max(templateLib.gapTrain)+cfg.rawGapSearchMargin);
z0 = min(max(z0,lb),ub);

penalty = 10*max(std(V),1e-3);
residual = @(z) per_sensor_residual(z,t,V,x,sensorId,templateLib,penalty);
startBank = z0;
if S == 3
    startBank = [startBank; z0 + [0.04, -0.02, -0.04, zeros(1,numel(z0)-3)]];
    startBank = [startBank; z0 + [-0.04, 0.02, 0.04, zeros(1,numel(z0)-3)]];
end
bestSse = Inf; z = z0; solveInfo = struct();
for istart = 1:size(startBank,1)
    zStart = min(max(startBank(istart,:),lb),ub);
    [zTry, infoTry] = solve_lsq_bounded(residual,zStart,lb,ub, ...
        get_field(cfg,'perSensorMaxIter',300),1e-10,1e-10);
    rTry = residual(zTry); sseTry = dot(rTry,rTry);
    if sseTry < bestSse
        bestSse = sseTry; z = zTry; solveInfo = infoTry;
    end
end

gHat = z(1:S);
dx = z(S+1);
p = z(S+2:end);
VFit = eval_per_sensor_voltage(gHat,dx,p,t,x,sensorId,templateLib);
rr = VFit-V; valid=isfinite(rr);
R = shared;
R.method = "per_sensor_gap_extension";
R.gap_mode = "per_sensor_absolute";
R.sensor_values = sensorValues;
R.sensor_gap_used = gHat(:);
R.sensor_gap_increment = gHat(:)-lowState.gHat;
R.g_used = mean(gHat);
R.dx_used = dx;
R.fit = struct('g',gHat,'dx',dx,'p',p,'f',[p(3),p(6)], ...
    'A',[p(1),p(4)],'phi',[p(2),p(5)],'theta',z, ...
    'VFit',VFit,'rmse',sqrt(mean(rr(valid).^2)), ...
    'solve_info',solveInfo);
R.f_id = [p(3),p(6)];
R.A_id = [p(1),p(4)];
R.phi_id = [p(2),p(5)];
R.VFit = VFit;
R.rmse = sqrt(mean(rr(valid).^2));
R.query_safe = queryInfo;
R.shared_seed = shared;
end

function r = per_sensor_residual(z,t,V,x,sensorId,lib,penalty)
S = max(sensorId);
g = z(1:S); dx=z(S+1); p=z(S+2:end);
y = eval_per_sensor_voltage(g,dx,p,t,x,sensorId,lib);
r = V-y; r(~isfinite(r))=penalty;
end

function y = eval_per_sensor_voltage(g,dx,p,t,x,sensorId,lib)
u = fit_u(p,t);
y = nan(size(x));
for s=1:numel(g)
    idx=sensorId==s;
    if any(idx), y(idx)=eval_gap_template(lib,g(s),x(idx)-dx-u(idx)); end
end
end

function value = get_field(s,name,defaultValue)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), value=s.(name); else, value=defaultValue; end
end
