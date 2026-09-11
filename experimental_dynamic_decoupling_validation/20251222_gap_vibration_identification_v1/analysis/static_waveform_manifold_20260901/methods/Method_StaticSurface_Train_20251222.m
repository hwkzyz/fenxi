function model = Method_StaticSurface_Train_20251222(trainY, trainGaps, g0)
nX=size(trainY,1); C=nan(nX,3); A=[ones(numel(trainGaps),1),1./trainGaps(:),log(trainGaps(:)./g0)];
for ix=1:nX, y=trainY(ix,:).'; ok=isfinite(y); if nnz(ok)>=3, C(ix,:)=A(ok,:)\y(ok); end, end
model.coeff=C; model.g0=g0; model.methodId='static_surface_basis_v1';
end
