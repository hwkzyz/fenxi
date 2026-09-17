function y = eval_gap_path_response(templateLib, g0, mu, tau, zeta, xq)
%EVAL_GAP_PATH_RESPONSE  Evaluate a static response surface along a path.

xq = xq(:);
xLib = zeta .* (xq - tau);
gPath = g0 + mu .* (xq - tau);
y = zeros(size(xq));
for i = 1:numel(xq)
    y(i) = eval_gap_template(templateLib, gPath(i), xLib(i));
end
end
