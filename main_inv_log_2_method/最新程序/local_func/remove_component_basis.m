function rPerp = remove_component_basis(r, B)
%remove_component_basis  Remove the component of r along the column space of B.

r = r(:);
B = B(:,:);
if isempty(B)
    rPerp = r;
    return;
end

keep = vecnorm(B, 2, 1) > 0;
B = B(:, keep);
if isempty(B)
    rPerp = r;
    return;
end

coef = B \ r;
rPerp = r - B * coef;
end
