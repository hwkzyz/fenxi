function [y, H] = build_varpro_design(V, x, t, templateLib, gHat, f1, f2)
%build_varpro_design  Build linearized variable-projection system.

F0 = eval_gap_template(templateLib, gHat, x);
Fx = eval_gap_derivative(templateLib, gHat, x);
s1 = sin(2*pi*f1*t);
c1 = cos(2*pi*f1*t);
s2 = sin(2*pi*f2*t);
c2 = cos(2*pi*f2*t);

y = V - F0;
H = [-Fx, -Fx .* s1, -Fx .* c1, -Fx .* s2, -Fx .* c2];
end
