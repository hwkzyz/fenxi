function info=compute_dx_vibration_identifiability(highMap,lib,g,dx,A,phi,f)
%COMPUTE_DX_VIBRATION_IDENTIFIABILITY  Local dx/AB column correlations.
t=highMap.t_v(:);x=highMap.x_v(:);A=A(:).';phi=phi(:).';f=f(:).';u=zeros(size(t));
for k=1:numel(f),u=u+A(k)*sin(2*pi*f(k)*t+phi(k));end
Fx=eval_gap_derivative(lib,g,x-dx-u);J=-Fx;names="dx";
for k=1:numel(f)
 J=[J,-Fx.*sin(2*pi*f(k)*t),-Fx.*cos(2*pi*f(k)*t)]; %#ok<AGROW>
 names=[names,"a"+k,"b"+k]; %#ok<AGROW>
end
ok=all(isfinite(J),2);J=J(ok,:);scale=vecnorm(J);scale(scale<=eps)=1;Jn=J./scale;G=Jn'*Jn;
corrDx=G(1,2:end);info=struct('column_names',names,'normalized_gram',G, ...
 'dx_dynamic_correlation',corrDx,'max_abs_dx_dynamic_correlation',max(abs(corrDx)), ...
 'normalized_gram_condition',cond(G),'weak_dx_identifiability',max(abs(corrDx))>0.9);
end
