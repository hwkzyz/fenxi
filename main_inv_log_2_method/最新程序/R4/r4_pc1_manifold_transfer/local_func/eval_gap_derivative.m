function y = eval_gap_derivative(templateLib, g, xq, sensorId)
%eval_gap_derivative  Evaluate x-derivative of interpolated template.

if isfield(templateLib, 'fixedPathIncrement') && templateLib.fixedPathIncrement
    xq = xq(:);
    h = 1e-4;
    % The cache is dense in both gap and position and is the same forward
    % surface used by candidate replay. Differentiating it avoids thousands
    % of scalar interp1 calls during iterative VP updates while preserving
    % all waveform samples.
    if nargin<4,sensorId=[];end
    y = (eval_fixed_path_template(templateLib, g, xq + h, sensorId) - ...
        eval_fixed_path_template(templateLib, g, xq - h, sensorId)) ./ (2*h);
    if isfield(templateLib,'lowDerivativeBySensor')
        % Replace only the PCHIP derivative of the measured low-speed term;
        % the calibrated gap-increment surface remains differentiated by a
        % centered finite difference.
        lowNum=(low_only(templateLib,xq+h,sensorId)-low_only(templateLib,xq-h,sensorId))./(2*h);
        lowStable=derivative_only(templateLib,xq,sensorId);
        y=y-lowNum+lowStable;
    end
    return;
end

[~, dGrid] = eval_gap_grid(templateLib, g);
y = interp1(templateLib.xGrid, dGrid(:), xq, 'pchip', 'extrap');
end

function y=low_only(lib,xq,sensorId)
y=select_interp(lib.xTemplate,lib.templateLow,lib.sensorIds,xq,sensorId);
end
function y=derivative_only(lib,xq,sensorId)
y=select_interp(lib.xTemplate,lib.lowDerivativeBySensor,lib.sensorIds,xq,sensorId);
end
function y=select_interp(x,T,ids,xq,sensorId)
if size(T,2)==1,y=interp1(x,T,xq,'pchip','extrap');return;end
if nargin<5||isempty(sensorId),sensorId=ids(1);end
if isscalar(sensorId),[ok,c]=ismember(sensorId,ids);if ~ok,error('Unknown sensor.');end;y=interp1(x,T(:,c),xq,'pchip','extrap');return;end
y=nan(size(xq));sensorId=sensorId(:);for sid=unique(sensorId).',[ok,c]=ismember(sid,ids);if ~ok,error('Unknown sensor.');end;q=sensorId==sid;y(q)=interp1(x,T(:,c),xq(q),'pchip','extrap');end
end
