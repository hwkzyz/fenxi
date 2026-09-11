function path = R5_Select_StatePath_SlowPrior(profileCells, lambda)
%R5_SELECT_STATEPATH_SLOWPRIOR Dynamic-programming state path from profile costs.
arguments, profileCells cell; lambda (1,1) double {mustBeNonnegative} = 0; end
nW=numel(profileCells); assert(nW>0,'R5:StatePath','No profile windows.');
states=unique([profileCells{1}.stateProfiles.stateValue]); states=states(isfinite(states)); nS=numel(states);
cost=inf(nW,nS);
for iw=1:nW
    P=profileCells{iw}.stateProfiles;
    for is=1:nS
        q=find([P.stateValue]==states(is) & strcmp({P.status},'profile_fitted'),1);
        if ~isempty(q), cost(iw,is)=P(q).J; end
    end
end
cost=cost./max(median(cost(isfinite(cost)),'omitnan'),eps);
J=inf(nW,nS); back=zeros(nW,nS); J(1,:)=cost(1,:);
for iw=2:nW
    for is=1:nS
        [prior,k]=min(J(iw-1,:)+lambda*(states(is)-states).^2);
        J(iw,is)=cost(iw,is)+prior; back(iw,is)=k;
    end
end
[total,k]=min(J(end,:)); idx=zeros(nW,1); idx(end)=k;
for iw=nW:-1:2, idx(iw-1)=back(iw,idx(iw)); end
path=struct('stateValue',states(idx).','stateIndex',idx.','totalCost',total, ...
    'lambda',lambda,'profileCost',cost,'status','slow_prior_path');
end
