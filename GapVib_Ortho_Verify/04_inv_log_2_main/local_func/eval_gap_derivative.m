function y = eval_gap_derivative(templateLib, g, xq)
%eval_gap_derivative  Evaluate x-derivative of interpolated template.

[~, dGrid] = eval_gap_grid(templateLib, g);
y = interp1(templateLib.xGrid, dGrid(:), xq, 'pchip', 'extrap');
end
