function eta = eta_gap(templateLib, g)
%eta_gap  Identifiability index after removing the translation direction.

x = templateLib.xGrid(:);
Fx = eval_gap_derivative(templateLib, g, x);
h = 1e-3;
Fg = (eval_gap_template(templateLib, g + h, x) - ...
    eval_gap_template(templateLib, g - h, x)) / (2*h);
FgPerp = Fg - Fx * (dot(Fx, Fg) / max(dot(Fx, Fx), eps));
eta = norm(FgPerp) / max(norm(Fg), eps);
end
