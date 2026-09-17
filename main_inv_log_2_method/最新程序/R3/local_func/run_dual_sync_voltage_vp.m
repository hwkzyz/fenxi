function R=run_dual_sync_voltage_vp(map,lib,cfg,staticState)
%RUN_DUAL_SYNC_VOLTAGE_VP Integer-order-pair dual voltage solver.
solverTic=tic;
t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);fRot=cfg.RPM_high/60;
if isfield(map,'S_v') && numel(map.S_v)==numel(x),sensorId=map.S_v(:);else,sensorId=[];end
if isfield(map,'theta_v')&&numel(map.theta_v)==numel(t)
    theta=map.theta_v(:);angleSource="measured_revolution_timing";
else
    theta=2*pi*fRot*t;angleSource="nominal_constant_speed";
end
range=get_field(cfg,'route30DualFrequencyRangeHz',[300 1500]);
eoGrid=ceil(range(1)/fRot):floor(range(2)/fRot);pairs=nchoosek(eoGrid,2);
% A known excitation-order pair is a valid physical prior for a targeted
% resonance diagnostic. The normal application path still enumerates pairs.
fixedPairs=get_field(cfg,'dualSyncCandidateEOPairs',[]);
if ~isempty(fixedPairs)
    fixedPairs=sort(unique(round(fixedPairs), 'rows'),2);
    if size(fixedPairs,2)~=2 || any(fixedPairs(:,1)>=fixedPairs(:,2)) || ...
            any(~ismember(fixedPairs(:),eoGrid))
        error('dualSyncCandidateEOPairs must contain distinct EO pairs within the search range.');
    end
    pairs=fixedPairs;
end
strategy=lower(string(get_field(cfg,'dualSyncSearchStrategy','global_budget')));
if strategy=="funnel_v1"
    R=run_dual_sync_voltage_vp_funnel(map,lib,cfg,staticState);
    return
end
S=sin(theta*eoGrid);C=cos(theta*eoGrid);
validStrategies=["global_budget","balanced_per_pair","coarse_top_m_multistart"];
if ~ismember(strategy,validStrategies)
    error('Unknown dualSyncSearchStrategy: %s.',strategy);
end
gHalf=get_field(cfg,'route30GapHalfWidthMm',.70);gN=get_field(cfg,'dualSyncGapCount',21);
if isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement
    gGrid=linspace(staticState.gHat-gHalf,staticState.gHat+gHalf,gN);
else
    lo=max(.05,min(lib.gapTrain)-cfg.rawGapSearchMargin);hi=max(lib.gapTrain)+cfg.rawGapSearchMargin;
    gGrid=linspace(max(lo,staticState.gHat-gHalf),min(hi,staticState.gHat+gHalf),gN);
end
keepPerGap=max(1,round(get_field(cfg,'dualSyncKeepPerGap',100)));
bank=repmat(empty_candidate(),0,1);
coarseTic=tic;
for g=gGrid
    F0=eval_gap_template(lib,g,x,sensorId);Fx=eval_gap_derivative(lib,g,x,sensorId);y=V-F0;
    valid=isfinite(Fx)&isfinite(y);q=-Fx(valid);yv=y(valid);Sv=S(valid,:);Cv=C(valid,:);
    rows=repmat(empty_candidate(),size(pairs,1),1);
    for ip=1:size(pairs,1)
        ia=pairs(ip,1)-eoGrid(1)+1;ib=pairs(ip,2)-eoGrid(1)+1;
        X=[q,q.*Sv(:,ia),q.*Cv(:,ia),q.*Sv(:,ib),q.*Cv(:,ib)];
        beta=X\yv;r=yv-X*beta;
        rows(ip)=pack_candidate(g,pairs(ip,:),beta,mean(r.^2),fRot);
    end
    [~,ord]=sort([rows.sse],'ascend');
    if strategy=="global_budget"
        keep=ord(1:min(keepPerGap,numel(ord)));
    else
        % Balanced and two-stage strategies must not lose an EO pair before
        % their explicit per-pair allocation stage.
        keep=ord;
    end
    bank=[bank;rows(keep)]; %#ok<AGROW>
end
coarseTime=toc(coarseTic);
[~,ord]=sort([bank.sse],'ascend');iterBudget=max(1,round(get_field(cfg,'dualSyncIteratedReplayCount',2100)));
iterCount=min(iterBudget,numel(ord));
replayTic=tic;
for ii=1:iterCount
    bank(ord(ii))=iterate_candidate(bank(ord(ii)),theta,x,V,lib,sensorId,...
        get_field(cfg,'dualSyncIteratedReplayMaxIter',2));
end
replayTime=toc(replayTic);
[pairEO,pairScore,pairRank]=rank_pairs(bank,pairs);
refineIdx=select_refine_indices(bank,pairEO,pairRank,strategy,cfg);
refineCount=numel(refineIdx);
best=[];bestRmse=Inf;refinedFrequency=zeros(0,2);refinedRmse=[];refinedFit=cell(0,1);
refineTic=tic;
for ii=1:refineCount
    c=bank(refineIdx(ii));z0=[c.g,c.dx,c.coef(:).'];
    gh=get_field(cfg,'route30SingleGapHalfWidthMm',.12);dh=get_field(cfg,'route30DxHalfWidthMm',.5);
    au=get_field(cfg,'route30AmplitudeUpperMm',1.5);
    lb=[max(.05,c.g-gh),-dh,-au*ones(1,4)];ub=[c.g+gh,dh,au*ones(1,4)];
    if ~(isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement)
        lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);
        ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
    end
    [z,si]=solve_lsq_bounded(@(v)fixed_pair_residual(v,c.eo,theta,V,x,lib,sensorId),...
        z0,lb,ub,get_field(cfg,'dualSyncMaxIter',180),1e-12,1e-12);
    [fit,VFit,rmse]=pack_fit(z,c.eo,theta,V,x,lib,sensorId,fRot,si);
    refinedFrequency(end+1,:)=fit.f;refinedRmse(end+1,1)=rmse; %#ok<AGROW>
    fit.VFit=VFit;refinedFit{end+1,1}=fit; %#ok<AGROW>
    if rmse<bestRmse,best=fit;best.VFit=VFit;bestRmse=rmse;end
end
refineTime=toc(refineTic);
amplitudeFloor=get_field(cfg,'structuredComponentAmplitudeFloorMm',.075);
eligible=cellfun(@(q)min(q.A)>=amplitudeFloor,refinedFit);
if any(eligible)
    idx=find(eligible);
    [~,j]=min(refinedRmse(idx));best=refinedFit{idx(j)};
end
R=struct('method',"dual_sync_voltage_vp",'model_order',2,'mode',"dual_sync_sync",...
    'g_used',best.g,'dx_used',best.dx,'A_id',best.A,'phi_id',best.phi,...
    'f_id',best.f,'eo_id',best.eo,'rmse',best.rmse,'VFit',best.VFit,...
    'fit',best,'angle_source',angleSource,'eo_grid',eoGrid(:),...
    'candidate_count',numel(bank),'used_sample_count',numel(t),...
    'all_trusted_samples_used',true,'staticState',staticState,...
    'component_amplitude_floor_mm',amplitudeFloor,...
    'search_strategy',strategy,'coarse_candidate_count',numel(bank),...
    'iterated_replay_count',iterCount,'nonlinear_refine_count',refineCount,...
    'coarse_time_s',coarseTime,'replay_time_s',replayTime,...
    'refine_time_s',refineTime,'solver_time_s',toc(solverTic));
C=frequency_solution_confidence(best.f,best.rmse,refinedFrequency(eligible,:),...
    refinedRmse(eligible),numel(t),cfg);
R=merge_struct(R,C);
if get_field(cfg,'returnDualSyncDiagnostics',false)
    [~,bo]=sort([bank.sse],'ascend');bs=bank(bo);
    R.diagnostic_bank_eo=vertcat(bs.eo);
    R.diagnostic_bank_sse=[bs.sse].';
    R.diagnostic_refined_eo=round(refinedFrequency/fRot);
    R.diagnostic_refined_rmse=refinedRmse;
    R.diagnostic_refined_eligible=eligible(:);
    R.diagnostic_pair_eo=pairEO;
    R.diagnostic_pair_score=pairScore;
    R.diagnostic_pair_rank=pairRank;
    R.diagnostic_refine_start_eo=vertcat(bank(refineIdx).eo);
    R.diagnostic_refine_start_sse=[bank(refineIdx).sse].';
    R.diagnostic_refine_start_g=[bank(refineIdx).g].';
    R.diagnostic_refine_start_dx=[bank(refineIdx).dx].';
end
end

function [pairEO,pairScore,pairRank]=rank_pairs(bank,pairs)
pairEO=pairs;pairScore=inf(size(pairs,1),1);bankEO=vertcat(bank.eo);score=[bank.sse].';
for ip=1:size(pairs,1)
    idx=all(bankEO==pairs(ip,:),2);
    if any(idx),pairScore(ip)=min(score(idx));end
end
[~,ord]=sort(pairScore,'ascend');pairRank=zeros(size(pairScore));pairRank(ord)=1:numel(ord);
end

function idx=select_refine_indices(bank,pairEO,pairRank,strategy,cfg)
score=[bank.sse].';bankEO=vertcat(bank.eo);
switch strategy
    case "global_budget"
        budget=max(1,round(get_field(cfg,'dualSyncRefineCount',50)));
        [~,ord]=sort(score,'ascend');idx=ord(1:min(budget,numel(ord)));
    case "balanced_per_pair"
        k=max(1,round(get_field(cfg,'dualSyncBalancedStartsPerPair',1)));
        idx=zeros(0,1);
        for ip=1:size(pairEO,1)
            q=find(all(bankEO==pairEO(ip,:),2));
            idx=[idx;stratified_indices(bank,q,score,k)]; %#ok<AGROW>
        end
    case "coarse_top_m_multistart"
        m=max(1,round(get_field(cfg,'dualSyncTopPairCount',10)));
        k=max(1,round(get_field(cfg,'dualSyncStartsPerTopPair',5)));
        chosen=find(pairRank<=min(m,numel(pairRank)));
        idx=zeros(0,1);
        for jj=1:numel(chosen)
            ip=chosen(jj);q=find(all(bankEO==pairEO(ip,:),2));
            idx=[idx;stratified_indices(bank,q,score,k)]; %#ok<AGROW>
        end
end
end

function idx=stratified_indices(bank,q,score,k)
% Retain distinct gap basins instead of spending K starts on adjacent gaps.
[~,go]=sort([bank(q).g],'ascend');q=q(go);n=numel(q);k=min(k,n);
% Use evenly spaced gap strata including both boundaries. The observed
% failure basin can sit at a gap boundary even when its coarse SSE is weak.
pos=unique(round(linspace(1,n,k)),'stable');idx=q(pos(:));
% If repeated rounded positions reduced the count, fill from the best
% remaining coarse candidates without duplicating an initial point.
if numel(idx)<k
    remain=setdiff(q,idx,'stable');[~,ord]=sort(score(remain),'ascend');
    idx=[idx;remain(ord(1:k-numel(idx)))];
end
end

function a=merge_struct(a,b)
names=fieldnames(b);for i=1:numel(names),a.(names{i})=b.(names{i});end
end

function c=empty_candidate()
c=struct('g',NaN,'dx',NaN,'eo',[NaN NaN],'f',[NaN NaN],...
    'coef',nan(4,1),'sse',Inf);
end

function c=pack_candidate(g,eo,beta,sse,fRot)
c=struct('g',g,'dx',beta(1),'eo',eo,'f',eo*fRot,'coef',beta(2:5),'sse',sse);
end

function c=iterate_candidate(c,theta,x,V,lib,sensorId,maxIter)
dx=c.dx;coef=c.coef(:);S1=sin(c.eo(1)*theta);C1=cos(c.eo(1)*theta);
S2=sin(c.eo(2)*theta);C2=cos(c.eo(2)*theta);
for it=1:maxIter
    u=coef(1)*S1+coef(2)*C1+coef(3)*S2+coef(4)*C2;
    z=x-dx-u;F0=eval_gap_template(lib,c.g,z,sensorId);Fx=eval_gap_derivative(lib,c.g,z,sensorId);
    valid=isfinite(V)&isfinite(F0)&isfinite(Fx);q=-Fx(valid);y=V(valid)-F0(valid);
    X=[q,q.*S1(valid),q.*C1(valid),q.*S2(valid),q.*C2(valid)];d=X\y;
    dx=dx+d(1);coef=coef+d(2:5);
end
u=coef(1)*S1+coef(2)*C1+coef(3)*S2+coef(4)*C2;
VFit=eval_gap_template(lib,c.g,x-dx-u,sensorId);r=V-VFit;good=isfinite(r);
c.dx=dx;c.coef=coef;c.sse=mean(r(good).^2);
end

function r=fixed_pair_residual(z,eo,theta,V,x,lib,sensorId)
u=z(3)*sin(eo(1)*theta)+z(4)*cos(eo(1)*theta)+...
  z(5)*sin(eo(2)*theta)+z(6)*cos(eo(2)*theta);
r=V-eval_gap_template(lib,z(1),x-z(2)-u,sensorId);r(~isfinite(r))=10*max(std(V),1e-3);
end

function [fit,VFit,rmse]=pack_fit(z,eo,theta,V,x,lib,sensorId,fRot,si)
A=[hypot(z(3),z(4)),hypot(z(5),z(6))];
phi=[atan2(z(4),z(3)),atan2(z(6),z(5))];
u=A(1)*sin(eo(1)*theta+phi(1))+A(2)*sin(eo(2)*theta+phi(2));
VFit=eval_gap_template(lib,z(1),x-z(2)-u,sensorId);r=V-VFit;good=isfinite(r);rmse=sqrt(mean(r(good).^2));
fit=struct('g',z(1),'dx',z(2),'A',A,'phi',phi,'eo',eo,'f',eo*fRot,...
    'theta',z,'rmse',rmse,'solve_info',si);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
