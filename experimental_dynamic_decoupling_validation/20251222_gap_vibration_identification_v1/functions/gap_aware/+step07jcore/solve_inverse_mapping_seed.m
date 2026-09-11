function [seedTable, info] = solve_inverse_mapping_seed(inv, bundle, cfg)
%SOLVE_INVERSE_MAPPING_SEED Fit integer-EO displacement candidates.

eoMin = ceil(cfg.freqSearchHz(1)/bundle.rotFreqMeanHz);
eoMax = floor(cfg.freqSearchHz(2)/bundle.rotFreqMeanHz);
eoList = max(1,eoMin):max(eoMin,eoMax);
u = inv.uInv(:); th = bundle.Theta(:); w = inv.weight(:);
valid = inv.valid(:) & isfinite(u) & isfinite(th) & isfinite(w) & w>0;
rows = repmat(struct('EO',NaN,'A',NaN,'phi',NaN,'dx',NaN,'weightedRmseMm',Inf,'validPointCount',0,'candidateRho',NaN), numel(eoList),1);
for i=1:numel(eoList)
    eo=eoList(i); H=[ones(nnz(valid),1),sin(eo*th(valid)),cos(eo*th(valid))];
    ww=w(valid); yy=u(valid); R=diag(sqrt(ww)); beta=(R*H)\(R*yy);
    res=yy-H*beta; rows(i).EO=eo; rows(i).dx=beta(1); rows(i).A=hypot(beta(2),beta(3)); rows(i).phi=atan2(beta(3),beta(2)); rows(i).weightedRmseMm=sqrt(sum(ww.*res.^2)/sum(ww)); rows(i).validPointCount=nnz(valid);
end
seedTable=struct2table(rows); seedTable=sortrows(seedTable,{'weightedRmseMm','EO'},{'ascend','ascend'});
if height(seedTable)>=2, rho=(seedTable.weightedRmseMm(2)-seedTable.weightedRmseMm(1))/max(seedTable.weightedRmseMm(2),eps); else, rho=NaN; end
info=struct('candidate_rho',rho,'candidate_count',height(seedTable),'eo_min',eoMin,'eo_max',eoMax,'valid_fraction',mean(valid));
end
