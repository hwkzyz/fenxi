function fit = pack_nonlinear_fit(dx, p, t)
%pack_nonlinear_fit  Normalize nonlinear fit parameters.

[fSort, order] = sort([p(3), p(6)]);
A = [abs(p(1)), abs(p(4))];
phi = [wrap_pi_local(p(2)), wrap_pi_local(p(5))];
fit.dx = dx;
fit.f = fSort;
fit.A = A(order);
fit.phi = phi(order);
fit.p = [fit.A(1), fit.phi(1), fit.f(1), fit.A(2), fit.phi(2), fit.f(2)];
fit.delta_sample = fit.dx + fit_u(fit.p, t);
end
