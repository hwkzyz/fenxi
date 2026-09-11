function [inv, summary] = inverse_map_experimental_displacement(bundle, dgRefMm, dxRefMm, model, cfg)
%INVERSE_MAP_EXPERIMENTAL_DISPLACEMENT Bounded two-branch inversion of F_s.
% model.sensor(is) must provide xGrid and evaluate(dgMm,xQuery), returning
% voltage in mV and, optionally, dVoltageDx in mV/mm. No extrapolation is
% permitted; the caller owns the exact experimental forward evaluator.

arguments
    bundle (1,1) struct
    dgRefMm (:,1) double
    dxRefMm (1,1) double
    model struct
    cfg (1,1) struct
end

V = local_field(bundle, {'V','v'});
X = local_field(bundle, {'X','x'});
W = local_field(bundle, {'W','fit_weight'}, true);
Theta = local_field(bundle, {'Theta','theta'}, true);
sensorIndex = local_field(bundle, {'sensorIndex','sensor_index'});
sensorIds = local_field(bundle, {'sensorIds','sensor_ids'});
V = V(:); X = X(:); sensorIndex = sensorIndex(:);
if isempty(W), W = ones(size(V)); else, W = W(:); end
if isempty(Theta), Theta = nan(size(V)); else, Theta = Theta(:); end
nSensor = numel(sensorIds);
if numel(dgRefMm) == 1, dgRefMm = repmat(dgRefMm,nSensor,1); end
if isfield(model,'sensor'), sensorModels = model.sensor; else, sensorModels = model; end
if numel(dgRefMm) ~= nSensor || numel(sensorModels) ~= nSensor
    error('inverse_map_experimental_displacement:SensorCount', ...
        'dgRefMm and model.sensor must match bundle sensor count.');
end

tolV = local_cfg(cfg, 'inverseVoltageToleranceMv', 1e-8);
tolX = local_cfg(cfg, 'inverseXToleranceMm', 1e-9);
minSlope = local_cfg(cfg, 'inverseMinDerivativeMvPerMm', 1e-6);
maxIter = max(20, round(local_cfg(cfg, 'inverseMaxBisectionIter', 80)));

inv = struct('xInv',nan(size(X)), 'uInv',nan(size(X)), ...
    'branchId',zeros(size(X)), 'valid',false(size(X)), ...
    'querySafe',false(size(X)), 'boundaryHit',false(size(X)), ...
    'noRoot',true(size(X)), 'voltageResidualMv',nan(size(X)), ...
    'dVoltageDx',nan(size(X)), 'weight',nan(size(X)), ...
    'sensorIds',sensorIds(:), 'dgRefMm',dgRefMm(:), 'dxRefMm',dxRefMm);

for is = 1:nSensor
    mask = sensorIndex == is;
    if ~any(mask), continue; end
    s = sensorModels(is);
    xGrid = double(s.xGrid(:));
    xGrid = xGrid(isfinite(xGrid));
    if numel(xGrid) < 3, continue; end
    xGrid = unique(sort(xGrid));
    [vGrid, dGrid] = local_eval(s, dgRefMm(is), xGrid);
    finiteGrid = isfinite(vGrid) & isfinite(xGrid);
    xGrid = xGrid(finiteGrid); vGrid = vGrid(finiteGrid);
    if numel(xGrid) < 3, continue; end
    [~, iPeak] = max(vGrid);
    [~, iValley] = min(vGrid);
    iExt = unique([iPeak iValley]);
    for k = find(mask(:)).'
        y = V(k);
        if ~isfinite(y), continue; end
        roots = [];
        branches = [];
        for ie = iExt(:).'
            ranges = {[1 ie], [ie numel(xGrid)]};
            for ib = 1:2
                lo = ranges{ib}(1); hi = ranges{ib}(2);
                if hi <= lo, continue; end
                xa = xGrid(lo:hi); ya = vGrid(lo:hi) - y;
                cross = find(ya(1:end-1).*ya(2:end) <= 0 & ...
                    isfinite(ya(1:end-1)) & isfinite(ya(2:end)), 1, 'first');
                if isempty(cross), continue; end
                a = xa(cross); b = xa(cross+1);
                fa = ya(cross); fb = ya(cross+1);
                if abs(fa) <= tolV, xr = a;
                elseif abs(fb) <= tolV, xr = b;
                else
                    for it = 1:maxIter
                        c = 0.5*(a+b); [vc,~] = local_eval(s,dgRefMm(is),c); fc = vc-y;
                        if abs(fc) <= tolV || abs(b-a) <= tolX, break; end
                        if fa*fc <= 0, b=c; fb=fc; else, a=c; fa=fc; end
                    end
                    xr = c;
                end
                if isempty(roots) || all(abs(roots-xr) > 10*tolX)
                    roots(end+1) = xr; %#ok<AGROW>
                    branches(end+1) = ib + 2*(ie ~= iPeak); %#ok<AGROW>
                end
            end
        end
        if isempty(roots), continue; end
        [~, pick] = min(abs((X(k)-dxRefMm) - roots));
        xr = roots(pick); [vr, dr] = local_eval(s,dgRefMm(is),xr);
        inv.xInv(k) = xr;
        inv.uInv(k) = X(k)-dxRefMm-xr;
        inv.branchId(k) = branches(pick);
        inv.voltageResidualMv(k) = vr-y;
        inv.dVoltageDx(k) = dr;
        inv.querySafe(k) = xr >= min(xGrid)-tolX && xr <= max(xGrid)+tolX;
        inv.boundaryHit(k) = abs(xr-xGrid(1)) <= tolX || abs(xr-xGrid(end)) <= tolX;
        inv.noRoot(k) = false;
        inv.valid(k) = inv.querySafe(k) && isfinite(dr) && abs(dr) >= minSlope;
    end
end

inv.weight = max(W, local_cfg(cfg,'weightFloor',0.05)) .* abs(inv.dVoltageDx).^2;
inv.weight(~inv.valid | ~isfinite(inv.weight)) = 0;
summary = struct();
summary.inverse_valid_fraction = mean(inv.valid);
summary.inverse_voltage_rmse = sqrt(mean(inv.voltageResidualMv(inv.valid).^2,'omitnan'));
summary.query_safe_fraction = mean(inv.querySafe);
summary.branch_switch_rate = local_branch_switch_rate(inv);
summary.no_root_fraction = mean(inv.noRoot);
summary.boundary_fraction = mean(inv.boundaryHit);
end

function [v,dv] = local_eval(s,dg,x)
[v, aux] = s.evaluate(dg,x);
v = double(v(:));
if isstruct(aux) && isfield(aux,'dVoltageDx'), dv = double(aux.dVoltageDx(:));
else, dv = nan(size(v)); end
if numel(dv) ~= numel(v), dv = nan(size(v)); end
end

function v = local_field(s,names,optional)
optional = nargin > 2 && optional; v = [];
for i=1:numel(names), if isfield(s,names{i}), v=s.(names{i}); return; end, end
if ~optional, error('inverse_map_experimental_displacement:Field','Missing bundle field.'); end
end

function v = local_cfg(cfg,name,default)
if isfield(cfg,name) && isfinite(cfg.(name)), v=cfg.(name); else, v=default; end
end

function r = local_branch_switch_rate(inv)
r = NaN; ids = inv.branchId(inv.valid); if numel(ids)>1, r=mean(ids(2:end)~=ids(1:end-1)); end
end
