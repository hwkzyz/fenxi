function y = eval_fixed_path_template(templateLib, gQuery, xq, sensorId)
% Evaluate the fixed low-speed path plus an absolute high-speed gap query.
pc = templateLib.pathCal;
querySize = size(xq);
xq = xq(:);
xLib = pc.zeta .* (xq - pc.tau);
gBase = pc.g0 + pc.mu .* (xq - pc.tau);
delta = gQuery - pc.g0;
gQueryPath = gBase + delta;
if isfield(templateLib, 'lowInterpolant') && ~isempty(templateLib.lowInterpolant)
    Tlow = templateLib.lowInterpolant(xq);
else
    Tlow = interpolate_sensor_template(templateLib,xq,sensorId);
end

if isfield(templateLib, 'pathCache') && ~isempty(templateLib.pathCache)
    C = templateLib.pathCache;
    if isfield(C, 'interpolant') && ~isempty(C.interpolant)
        yBase = C.interpolant(gBase, xLib);
        yQuery = C.interpolant(gQueryPath, xLib);
    else
        yBase = interp2(C.xGrid, C.gGrid, C.SGrid, xLib, gBase, 'linear');
        yQuery = interp2(C.xGrid, C.gGrid, C.SGrid, xLib, gQueryPath, 'linear');
    end
else
    if max(gBase)-min(gBase) < 1e-12
        yBase = eval_gap_template(templateLib.baseLib, gBase(1), xLib);
    else
        yBase = eval_variable_base(templateLib.baseLib, gBase, xLib);
    end
    if max(gQueryPath)-min(gQueryPath) < 1e-12
        yQuery = eval_gap_template(templateLib.baseLib, gQueryPath(1), xLib);
    else
        yQuery = eval_variable_base(templateLib.baseLib, gQueryPath, xLib);
    end
end
y = reshape(Tlow + pc.kappa .* (yQuery - yBase),querySize);
end

function Tlow=interpolate_sensor_template(templateLib,xq,sensorId)
T=templateLib.templateLow;
if size(T,2)==1
    Tlow=interp1(templateLib.xTemplate,T,xq,'pchip','extrap');
    return;
end
if nargin<3 || isempty(sensorId),sensorId=templateLib.sensorIds(1);end
if isscalar(sensorId)
    [hit,col]=ismember(sensorId,templateLib.sensorIds);
    if ~hit,error('No low-speed template is available for sensor %g.',sensorId);end
    Tlow=interp1(templateLib.xTemplate,T(:,col),xq,'pchip','extrap');
    return;
end
sensorId=sensorId(:);
if numel(sensorId)~=numel(xq),error('sensorId must be scalar or match xq.');end
Tlow=nan(size(xq));
for sid=unique(sensorId(:)).'
    [hit,col]=ismember(sid,templateLib.sensorIds);
    if ~hit,error('No low-speed template is available for sensor %g.',sid);end
    idx=sensorId==sid;
    Tlow(idx)=interp1(templateLib.xTemplate,T(:,col),xq(idx),'pchip','extrap');
end
end

function y = eval_variable_base(baseLib, g, x)
g = g(:); x = x(:); y = nan(size(g));
for i = 1:numel(g)
    y(i) = eval_gap_template(baseLib, g(i), x(i));
end
end
