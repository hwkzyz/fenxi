function precomp = build_wave_varpro_precompute(t, V, x, templateLib, gHat, dxHat, cfg)
%build_wave_varpro_precompute  Precompute reusable terms for fixed-gap VARPRO.

t = t(:);
V = V(:);
x = x(:);

F0 = eval_gap_template(templateLib, gHat, x - dxHat);
Fx = eval_gap_derivative(templateLib, gHat, x - dxHat);
DV = V - F0;

precomp.F0 = F0;
precomp.Fx = Fx;
precomp.DV = DV;
precomp.S1 = sin(2*pi*t*cfg.f1Grid(:).');
precomp.C1 = cos(2*pi*t*cfg.f1Grid(:).');
precomp.S2 = sin(2*pi*t*cfg.f2Grid(:).');
precomp.C2 = cos(2*pi*t*cfg.f2Grid(:).');
end
