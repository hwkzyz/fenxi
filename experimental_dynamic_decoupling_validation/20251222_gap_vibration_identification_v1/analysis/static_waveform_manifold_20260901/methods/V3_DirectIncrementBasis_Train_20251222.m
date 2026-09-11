function model = V3_DirectIncrementBasis_Train_20251222(trainY, gaps)
nX=size(trainY,1); nB=size(trainY,2); nG=numel(gaps); A=[];
for il=1:nG-1, it=il+1; A(end+1,:)=[1,1/gaps(it)-1/gaps(il),log(gaps(it)/gaps(il))]; end
AA=repmat(A,nB,1); C=nan(nX,3);
for ix=1:nX, y=[]; for il=1:nG-1, y=[y; squeeze(trainY(ix,:,il+1)-trainY(ix,:,il)).']; end; ok=isfinite(y); if nnz(ok)>=3, C(ix,:)=AA(ok,:)\y(ok); end, end
model.coeff=C; model.methodId='direct_increment_basis_v1';
end
