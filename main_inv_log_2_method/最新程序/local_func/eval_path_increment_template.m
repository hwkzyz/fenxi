function y = eval_path_increment_template(pathModel, deltaGap, xq, sensorId)
%EVAL_PATH_INCREMENT_TEMPLATE  Evaluate the fixed-path increment model.
%
% The low-speed measured template is retained as the baseline.  The static
% library contributes only the response increment caused by deltaGap:
%   T_low(x) + kappa*(R(g_base(x)+deltaGap,x_lib)-R(g_base(x),x_lib)).

if ~isscalar(deltaGap) || ~isfinite(deltaGap)
    error('deltaGap must be a finite scalar.');
end
xq = xq(:);
if nargin<4 || isempty(sensorId),sensorId=pathModel.sensorIds(1);end
pc = pathModel.pathCal;
xLocal = xq;
xLib = pc.zeta .* (xLocal - pc.tau);
gBase = pc.g0 + pc.mu .* (xLocal - pc.tau);
gQuery = gBase + deltaGap;

rBase = eval_variable_gap_template(pathModel.templateLib, gBase, xLib);
rQuery = eval_variable_gap_template(pathModel.templateLib, gQuery, xLib);
tLow = interp1(pathModel.xTemplate, select_template(pathModel,sensorId), xLocal, ...
    'pchip', 'extrap');
y = tLow + pc.kappa .* (rQuery - rBase);
end

function t=select_template(pathModel,sensorId)
if size(pathModel.templateLow,2)==1,t=pathModel.templateLow;return;end
[hit,col]=ismember(sensorId,pathModel.sensorIds);
if ~hit,error('No low-speed template is available for sensor %g.',sensorId);end
t=pathModel.templateLow(:,col);
end

function y = eval_variable_gap_template(templateLib, g, xq)
g = g(:);
xq = xq(:);
if numel(g) ~= numel(xq)
    error('Variable gap and query coordinate must have equal lengths.');
end
y = zeros(size(g));
for i = 1:numel(g)
    y(i) = eval_gap_template(templateLib, g(i), xq(i));
end
end
