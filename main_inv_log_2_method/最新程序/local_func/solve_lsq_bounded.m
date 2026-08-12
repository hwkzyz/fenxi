function [theta, info] = solve_lsq_bounded(resFun, theta0, lb, ub, maxIter, funcTol, stepTol)
%solve_lsq_bounded  Bound-constrained least-squares with fminsearch fallback.

theta0 = min(max(theta0, lb), ub);
if nargin < 6 || isempty(funcTol)
    funcTol = 1e-10;
end
if nargin < 7 || isempty(stepTol)
    stepTol = 1e-10;
end
info = struct('solver', "", 'iterations', NaN, 'funcCount', NaN, ...
    'exitflag', NaN, 'message', "");
if exist('lsqnonlin', 'file') == 2
    opts = optimoptions('lsqnonlin', 'Display', 'none', ...
        'FunctionTolerance', funcTol, 'StepTolerance', stepTol, ...
        'MaxIterations', maxIter, ...
        'MaxFunctionEvaluations', max(2000, 50 * maxIter));
    [theta, ~, ~, exitflag, output] = lsqnonlin(resFun, theta0, lb, ub, opts);
    info.solver = "lsqnonlin";
    info.iterations = get_output_field(output, 'iterations');
    info.funcCount = get_output_field(output, 'funcCount');
    info.exitflag = exitflag;
    info.message = string(get_output_field(output, 'message'));
else
    width = max(ub - lb, eps);
    z0 = log((theta0 - lb + 1e-9) ./ max(ub - theta0 + 1e-9, 1e-9));
    obj = @(z) sum(resFun(lb + width ./ (1 + exp(-z))).^2);
    [z, ~, exitflag, output] = fminsearch(obj, z0, optimset('Display', 'off', 'MaxIter', 1200));
    theta = lb + width ./ (1 + exp(-z));
    info.solver = "fminsearch";
    info.iterations = get_output_field(output, 'iterations');
    info.funcCount = get_output_field(output, 'funcCount');
    info.exitflag = exitflag;
    info.message = string(get_output_field(output, 'message'));
end
end

function val = get_output_field(S, name)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    val = S.(name);
else
    val = NaN;
end
end
