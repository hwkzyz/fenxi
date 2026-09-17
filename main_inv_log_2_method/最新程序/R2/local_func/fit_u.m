function u = fit_u(p, t)
%fit_u  Two-frequency displacement model.

u = p(1) * sin(2*pi*p(3)*t + p(2)) + ...
    p(4) * sin(2*pi*p(6)*t + p(5));
end
