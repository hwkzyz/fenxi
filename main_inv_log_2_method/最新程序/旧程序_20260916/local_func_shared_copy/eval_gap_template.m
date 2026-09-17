function y = eval_gap_template(templateLib, g, xq, sensorId)
%eval_gap_template  Evaluate interpolated static-gap template.

if isfield(templateLib, 'fixedPathIncrement') && templateLib.fixedPathIncrement
    if nargin<4,sensorId=[];end
    y = eval_fixed_path_template(templateLib, g, xq, sensorId);
    return;
end

[curveGrid, ~] = eval_gap_grid(templateLib, g);
y = interp1(templateLib.xGrid, curveGrid(:), xq, 'pchip', 'extrap');
end
