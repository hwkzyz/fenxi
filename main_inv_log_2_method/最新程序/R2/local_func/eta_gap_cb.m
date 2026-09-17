function eta = eta_gap_cb(templateLib, gap, modeName)
%eta_gap_cb  Gap identifiability after removing nuisance directions.

x = templateLib.xGrid(:);
F = eval_gap_template(templateLib, gap, x);
Fx = eval_gap_derivative(templateLib, gap, x);
h = 1e-3;
Fg = (eval_gap_template(templateLib, gap + h, x) - ...
    eval_gap_template(templateLib, gap - h, x)) / (2*h);

switch string(modeName)
    case "main"
        B = Fx;
    case "b"
        B = [ones(size(x)), Fx];
    case "cb"
        B = [ones(size(x)), F, Fx];
    otherwise
        error('Unknown modeName: %s', modeName);
end

FgPerp = remove_component_basis(Fg, B);
eta = norm(FgPerp) / max(norm(Fg), eps);
end
