function y = eval_gap_template(templateLib, g, xq)
%eval_gap_template  Evaluate interpolated static-gap template.

[curveGrid, ~] = eval_gap_grid(templateLib, g);
y = interp1(templateLib.xGrid, curveGrid(:), xq, 'pchip', 'extrap');
end
