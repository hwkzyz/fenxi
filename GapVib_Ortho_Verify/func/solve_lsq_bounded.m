function theta = solve_lsq_bounded(resFun, theta0, lb, ub, maxIter, funcTol, stepTol)
%solve_lsq_bounded  Bound-constrained least-squares with fminsearch fallback.

theta0 = min(max(theta0, lb), ub);
if nargin < 6 || isempty(funcTol)
    funcTol = 1e-10;
end
if nargin < 7 || isempty(stepTol)
    stepTol = 1e-10;
end
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', funcTol, 'StepTolerance', stepTol, ...
        'MaxIterations', maxIter);
    theta = lsqnonlin(resFun, theta0, lb, ub, opts);
else
    width = max(ub - lb, eps);
    z0 = log((theta0 - lb + 1e-9) ./ max(ub - theta0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    z = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 1200));
    theta = lb + width ./ (1 + exp(-z));
end
end
