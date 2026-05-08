function rPerp = remove_translation_component(r, Fx, lambda)
%remove_translation_component  Suppress residual component along translation.

if nargin < 3 || isempty(lambda)
    lambda = 1;
end
lambda = max(0, min(1, lambda));
rPerp = r - lambda * Fx * (dot(Fx, r) / max(dot(Fx, Fx), eps));
end
