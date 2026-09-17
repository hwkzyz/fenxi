function fit = pack_linear_fit(beta, freqPair, t)
%pack_linear_fit  Convert linear sine/cosine coefficients to A-phi form.

A = [hypot(beta(2), beta(3)), hypot(beta(4), beta(5))];
phi = [atan2(beta(3), beta(2)), atan2(beta(5), beta(4))];
[fSort, order] = sort(freqPair);
fit.dx = beta(1);
fit.f = fSort;
fit.A = A(order);
fit.phi = phi(order);
fit.p = [fit.A(1), fit.phi(1), fit.f(1), fit.A(2), fit.phi(2), fit.f(2)];
fit.delta_sample = fit.dx + fit_u(fit.p, t);
end
