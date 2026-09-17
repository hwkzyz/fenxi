function R=run_dual_sync_voltage_vp_funnel(map,lib,cfg,staticState)
%RUN_DUAL_SYNC_VOLTAGE_VP_FUNNEL Conservative dual-sync Funnel V1.
% All EO pairs survive the linearized VP stage.  Each pair contributes two
% separated gap states to exact replay and one fixed-gap damped GN step;
% twelve distinct EO pairs receive two joint short-GN steps and three
% distinct pairs receive the unchanged full bounded nonlinear refinement.
solverTic=tic;
t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);fRot=cfg.RPM_high/60;
assert_certified_map(x,cfg);
if isfield(map,'S_v')&&numel(map.S_v)==numel(x),sensorId=map.S_v(:);else,sensorId=[];end
if isfield(map,'theta_v')&&numel(map.theta_v)==numel(t)
    theta=map.theta_v(:);angleSource="measured_revolution_timing";
else
    theta=2*pi*fRot*t;angleSource="nominal_constant_speed";
end
range=get_field(cfg,'route30DualFrequencyRangeHz',[300 1500]);
eoGrid=ceil(range(1)/fRot):floor(range(2)/fRot);pairs=nchoosek(eoGrid,2);
fixedPairs=get_field(cfg,'dualSyncCandidateEOPairs',[]);
if ~isempty(fixedPairs)
    fixedPairs=sort(unique(round(fixedPairs),'rows'),2);
    if size(fixedPairs,2)~=2||any(fixedPairs(:,1)>=fixedPairs(:,2))||...
            any(~ismember(fixedPairs(:),eoGrid))
        error('dualSyncCandidateEOPairs must contain distinct EO pairs within the search range.');
    end
    pairs=fixedPairs;
end
S=sin(theta*eoGrid);C=cos(theta*eoGrid);
[gGrid,gBounds]=make_gap_grid(lib,cfg,staticState);

% Level 0: complete EO x gap linearized VP profile.
coarseTic=tic;bank=repmat(empty_candidate(),numel(gGrid)*size(pairs,1),1);ibank=0;
for ig=1:numel(gGrid)
    g=gGrid(ig);F0=eval_gap_template(lib,g,x,sensorId);
    Fx=eval_gap_derivative(lib,g,x,sensorId);y=V-F0;
    valid=isfinite(Fx)&isfinite(y);q=-Fx(valid);yv=y(valid);Sv=S(valid,:);Cv=C(valid,:);
    for ip=1:size(pairs,1)
        ia=pairs(ip,1)-eoGrid(1)+1;ib=pairs(ip,2)-eoGrid(1)+1;
        X=[q,q.*Sv(:,ia),q.*Cv(:,ia),q.*Sv(:,ib),q.*Cv(:,ib)];
        beta=X\yv;r=yv-X*beta;ibank=ibank+1;
        bank(ibank)=pack_candidate(g,pairs(ip,:),beta,mean(r.^2),fRot,ip,ig);
    end
end
coarseTime=toc(coarseTic);

% Level 1: two separated gap states for every EO pair, exact0 + fixed-g GN1.
level1Idx=zeros(2*size(pairs,1),1);k=0;
for ip=1:size(pairs,1)
    q=find([bank.pair_index].'==ip);score=[bank(q).sse_vp].';
    [~,j1]=min(score);blocked=abs([bank(q).gap_index].'-bank(q(j1)).gap_index)<=1;
    score2=score;score2(blocked)=Inf;
    if all(~isfinite(score2))
        [~,j2]=max(abs([bank(q).gap_index].'-bank(q(j1)).gap_index));
    else
        [~,j2]=min(score2);
    end
    k=k+1;level1Idx(k)=q(j1);k=k+1;level1Idx(k)=q(j2);
end
level1Tic=tic;
for i=1:numel(level1Idx)
    c=bank(level1Idx(i));c.sse_exact0=exact_sse(c,theta,x,V,lib,sensorId,cfg);
    [c,c.gn1_accepted]=fixed_gap_lm_step(c,theta,x,V,lib,sensorId,cfg);
    c.sse_gn1=exact_sse(c,theta,x,V,lib,sensorId,cfg);
    bank(level1Idx(i))=c;
end
level1Time=toc(level1Tic);

% Merge the two gap states and keep twelve distinct EO pairs.  Ten are
% selected by final GN1 residual; two optional rescue pairs require a
% top-half GN1 score and are ranked by relative improvement.
pairBest=zeros(size(pairs,1),1);pairScore=inf(size(pairs,1),1);
pairExact0=inf(size(pairs,1),1);pairDrop=-inf(size(pairs,1),1);
for ip=1:size(pairs,1)
    q=level1Idx([bank(level1Idx).pair_index].'==ip);
    s=[bank(q).sse_gn1];[pairScore(ip),j]=min(s);pairBest(ip)=q(j);
    pairExact0(ip)=bank(q(j)).sse_exact0;
    pairDrop(ip)=(pairExact0(ip)-pairScore(ip))/max(pairExact0(ip),eps);
end
[~,pairOrder]=sort(pairScore,'ascend');
useAllJointV2=get_field(cfg,'funnelAllPairsFirstJointStep',false);
joint1Idx=[];joint1Score=[];
level2Tic=tic;
joint1Time=0;joint2Time=0;
if useAllJointV2
    % V2: every discrete EO pair receives one gap-inclusive correction before pruning.
    joint1Idx=pairBest;
    joint1Tic=tic;
    for i=1:numel(joint1Idx)
        c=bank(joint1Idx(i));[c,ok]=joint_lm_step(c,theta,x,V,lib,sensorId,cfg,gBounds);
        c.short_gn_accepted=ok;c.sse_short=exact_sse(c,theta,x,V,lib,sensorId,cfg);bank(joint1Idx(i))=c;
    end
    joint1Time=toc(joint1Tic);
    joint1Score=[bank(joint1Idx).sse_short].';[~,joint1Order]=sort(joint1Score,'ascend');
    target=min(get_field(cfg,'funnelSecondJointPairCount',24),numel(joint1Order));
    shortIdx=joint1Idx(joint1Order(1:target));
    joint2Tic=tic;
    for i=1:numel(shortIdx)
        c=bank(shortIdx(i));[c,ok]=joint_lm_step(c,theta,x,V,lib,sensorId,cfg,gBounds);
        c.short_gn_accepted=c.short_gn_accepted+ok;c.sse_short=exact_sse(c,theta,x,V,lib,sensorId,cfg);bank(shortIdx(i))=c;
    end
    joint2Time=toc(joint2Tic);
else
    target=min(get_field(cfg,'funnelShortGnPairCount',40),numel(pairOrder));
    topDefault=max(1,target-2);topN=min(get_field(cfg,'funnelTopResidualPairCount',topDefault),target);
    chosen=pairOrder(1:topN);
    eligible=pairOrder(1:max(1,ceil(numel(pairOrder)/2)));eligible=setdiff(eligible,chosen,'stable');
    [~,dropOrder]=sort(pairDrop(eligible),'descend');rescue=eligible(dropOrder(1:min(target-numel(chosen),numel(dropOrder))));
    chosen=[chosen(:);rescue(:)];
    if numel(chosen)<target
        chosen=[chosen;setdiff(pairOrder,chosen,'stable')];chosen=chosen(1:target);
    end
    shortIdx=pairBest(chosen);shortSteps=get_field(cfg,'funnelJointGnSteps',2);joint2Tic=tic;
    for i=1:numel(shortIdx)
        c=bank(shortIdx(i));accepted=0;
        for it=1:shortSteps,[c,ok]=joint_lm_step(c,theta,x,V,lib,sensorId,cfg,gBounds);accepted=accepted+ok;end
        c.short_gn_accepted=accepted;c.sse_short=exact_sse(c,theta,x,V,lib,sensorId,cfg);bank(shortIdx(i))=c;
    end
    joint2Time=toc(joint2Tic);
end
level2Time=toc(level2Tic);
[~,shortOrder]=sort([bank(shortIdx).sse_short],'ascend');
fullCount=min(get_field(cfg,'funnelFullRefineCount',3),numel(shortOrder));
fullIdx=shortIdx(shortOrder(1:fullCount));

% Level 3: unchanged full bounded nonlinear solver, fixed at three starts
% during V1 validation so that the first/second/third rescue rates remain observable.
refineTic=tic;refinedFit=cell(fullCount,1);refinedRmse=inf(fullCount,1);
fullIterations=zeros(fullCount,1);fullFuncCount=zeros(fullCount,1);
for i=1:fullCount
    c=bank(fullIdx(i));z0=[c.g,c.dx,c.coef(:).'];
    gh=get_field(cfg,'route30SingleGapHalfWidthMm',.12);dh=get_field(cfg,'route30DxHalfWidthMm',.5);
    au=get_field(cfg,'route30AmplitudeUpperMm',1.5);
    lb=[max(gBounds(1),c.g-gh),-dh,-au*ones(1,4)];
    ub=[min(gBounds(2),c.g+gh), dh, au*ones(1,4)];
    [z,si]=solve_lsq_bounded(@(v)fixed_pair_residual(v,c.eo,theta,V,x,lib,sensorId,cfg),...
        z0,lb,ub,get_field(cfg,'dualSyncMaxIter',180),1e-12,1e-12);
    [fit,VFit,rmse]=pack_fit(z,c.eo,theta,V,x,lib,sensorId,fRot,si);
    fit.VFit=VFit;refinedFit{i}=fit;refinedRmse(i)=rmse;
    fullIterations(i)=si.iterations;fullFuncCount(i)=si.funcCount;
end
refineTime=toc(refineTic);
amplitudeFloor=get_field(cfg,'structuredComponentAmplitudeFloorMm',.075);
eligibleFit=cellfun(@(q)min(q.A)>=amplitudeFloor,refinedFit);
if any(eligibleFit),idx=find(eligibleFit);else,idx=(1:fullCount).';end
[~,j]=min(refinedRmse(idx));best=refinedFit{idx(j)};

methodName="dual_sync_voltage_vp_funnel_v1";strategyName="funnel_v1";
if useAllJointV2,methodName="dual_sync_voltage_vp_funnel_v2_all_joint1";strategyName="funnel_v2_all_joint1";end
R=struct('method',methodName,'model_order',2,...
    'mode',"dual_sync_sync",'g_used',best.g,'dx_used',best.dx,...
    'A_id',best.A,'phi_id',best.phi,'f_id',best.f,'eo_id',best.eo,...
    'rmse',best.rmse,'VFit',best.VFit,'fit',best,'angle_source',angleSource,...
    'eo_grid',eoGrid(:),'candidate_count',numel(bank),'used_sample_count',numel(t),...
    'all_trusted_samples_used',true,'staticState',staticState,...
    'component_amplitude_floor_mm',amplitudeFloor,'search_strategy',strategyName,...
    'coarse_candidate_count',numel(bank),'level1_candidate_count',numel(level1Idx),...
    'short_gn_pair_count',numel(shortIdx),'joint1_pair_count',numel(joint1Idx),...
    'nonlinear_refine_count',fullCount,...
    'coarse_time_s',coarseTime,'exact_gn1_time_s',level1Time,...
    'short_gn_time_s',level2Time,'joint1_time_s',joint1Time,'joint2_time_s',joint2Time,...
    'refine_time_s',refineTime,...
    'solver_time_s',toc(solverTic),'full_refine_iterations',fullIterations,...
    'full_refine_function_evaluations',fullFuncCount,...
    'level1_gn_accepted_count',sum([bank(level1Idx).gn1_accepted]),...
    'short_gn_accepted_count',sum([bank(shortIdx).short_gn_accepted]));
freq=zeros(fullCount,2);for i=1:fullCount,freq(i,:)=refinedFit{i}.f;end
Cconf=frequency_solution_confidence(best.f,best.rmse,freq,refinedRmse,numel(t),cfg);
R=merge_struct(R,Cconf);
if get_field(cfg,'returnDualSyncDiagnostics',false)
    R.diagnostic_pair_eo=pairs;R.diagnostic_pair_gn1_score=pairScore;
    R.diagnostic_pair_exact0_score=pairExact0;R.diagnostic_pair_relative_drop=pairDrop;
    if useAllJointV2
        R.diagnostic_joint1_eo=vertcat(bank(joint1Idx).eo);R.diagnostic_joint1_score=joint1Score;
    else
        R.diagnostic_joint1_eo=zeros(0,2);R.diagnostic_joint1_score=zeros(0,1);
    end
    R.diagnostic_level1_eo=vertcat(bank(level1Idx).eo);
    R.diagnostic_level1_gap=[bank(level1Idx).g].';
    R.diagnostic_level1_exact0=[bank(level1Idx).sse_exact0].';
    R.diagnostic_level1_gn1=[bank(level1Idx).sse_gn1].';
    R.diagnostic_short_eo=vertcat(bank(shortIdx).eo);
    R.diagnostic_short_sse=[bank(shortIdx).sse_short].';
    R.diagnostic_full_start_eo=vertcat(bank(fullIdx).eo);
    R.diagnostic_refined_eo=round(freq/fRot);R.diagnostic_refined_rmse=refinedRmse;
end
end

function [gGrid,bounds]=make_gap_grid(lib,cfg,state)
half=get_field(cfg,'route30GapHalfWidthMm',.70);n=get_field(cfg,'dualSyncGapCount',11);
if isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement
    bounds=[max(.05,state.gHat-half),state.gHat+half];
else
    bounds=[max(.05,min(lib.gapTrain)-cfg.rawGapSearchMargin),...
        max(lib.gapTrain)+cfg.rawGapSearchMargin];
    bounds=[max(bounds(1),state.gHat-half),min(bounds(2),state.gHat+half)];
end
gGrid=linspace(bounds(1),bounds(2),n);
end

function c=empty_candidate()
c=struct('g',NaN,'dx',NaN,'eo',[NaN NaN],'f',[NaN NaN],...
    'coef',nan(4,1),'sse_vp',Inf,'sse_exact0',Inf,'sse_gn1',Inf,...
    'sse_short',Inf,'pair_index',0,'gap_index',0,'gn1_accepted',false,...
    'short_gn_accepted',0);
end

function c=pack_candidate(g,eo,beta,sse,fRot,ip,ig)
c=empty_candidate();c.g=g;c.dx=beta(1);c.eo=eo;c.f=eo*fRot;
c.coef=beta(2:5);c.sse_vp=sse;c.pair_index=ip;c.gap_index=ig;
end

function s=exact_sse(c,theta,x,V,lib,sensorId,cfg)
u=harmonic_displacement(c,theta);
if ~candidate_admissible(c,x,u,cfg),s=Inf;return;end
fit=eval_gap_template(lib,c.g,x-c.dx-u,sensorId);
r=V-fit;if all(isfinite(r)),s=mean(r.^2);else,s=Inf;end
end

function u=harmonic_displacement(c,theta)
u=c.coef(1)*sin(c.eo(1)*theta)+c.coef(2)*cos(c.eo(1)*theta)+...
    c.coef(3)*sin(c.eo(2)*theta)+c.coef(4)*cos(c.eo(2)*theta);
end

function [c,accepted]=fixed_gap_lm_step(c,theta,x,V,lib,sensorId,cfg)
s1=sin(c.eo(1)*theta);c1=cos(c.eo(1)*theta);
s2=sin(c.eo(2)*theta);c2=cos(c.eo(2)*theta);
u=harmonic_displacement(c,theta);z=x-c.dx-u;
fit=eval_gap_template(lib,c.g,z,sensorId);fx=eval_gap_derivative(lib,c.g,z,sensorId);
r=V-fit;valid=isfinite(r)&isfinite(fx);J=[fx(valid),fx(valid).*s1(valid),...
    fx(valid).*c1(valid),fx(valid).*s2(valid),fx(valid).*c2(valid)];
p=[c.dx;c.coef(:)];sc=get_field(cfg,'funnelConditionalScaleMm',[.1 .1 .1 .1 .1]).';
lim=get_field(cfg,'funnelConditionalStepLimitMm',[.08 .12 .12 .12 .12]).';
[pNew,accepted]=damped_step(p,r(valid),J,sc,lim,@objective);
if accepted,c.dx=pNew(1);c.coef=pNew(2:5);end
    function val=objective(q)
        cc=c;cc.dx=q(1);cc.coef=q(2:5);val=exact_sse(cc,theta,x,V,lib,sensorId,cfg);
    end
end

function [c,accepted]=joint_lm_step(c,theta,x,V,lib,sensorId,cfg,gBounds)
s1=sin(c.eo(1)*theta);c1=cos(c.eo(1)*theta);
s2=sin(c.eo(2)*theta);c2=cos(c.eo(2)*theta);
u=harmonic_displacement(c,theta);z=x-c.dx-u;
fit=eval_gap_template(lib,c.g,z,sensorId);fx=eval_gap_derivative(lib,c.g,z,sensorId);
h=get_field(cfg,'funnelGapDerivativeStepMm',.002);
gm=max(gBounds(1),c.g-h);gp=min(gBounds(2),c.g+h);
fg=(eval_gap_template(lib,gp,z,sensorId)-eval_gap_template(lib,gm,z,sensorId))/max(gp-gm,eps);
r=V-fit;valid=isfinite(r)&isfinite(fx)&isfinite(fg);
J=[-fg(valid),fx(valid),fx(valid).*s1(valid),fx(valid).*c1(valid),...
    fx(valid).*s2(valid),fx(valid).*c2(valid)];
p=[c.g;c.dx;c.coef(:)];sc=get_field(cfg,'funnelJointScaleMm',[.08 .1 .1 .1 .1 .1]).';
lim=get_field(cfg,'funnelJointStepLimitMm',[.04 .08 .12 .12 .12 .12]).';
[pNew,accepted]=damped_step(p,r(valid),J,sc,lim,@objective);
if accepted,c.g=min(max(pNew(1),gBounds(1)),gBounds(2));c.dx=pNew(2);c.coef=pNew(3:6);end
    function val=objective(q)
        cc=c;cc.g=min(max(q(1),gBounds(1)),gBounds(2));cc.dx=q(2);cc.coef=q(3:6);
        val=exact_sse(cc,theta,x,V,lib,sensorId,cfg);
    end
end

function [pNew,accepted]=damped_step(p,r,J,scale,limit,objective)
old=objective(p);accepted=false;pNew=p;Js=J.*scale.';lambda=1e-3;
for trial=1:7
    A=Js.'*Js+lambda*diag(max(sum(Js.^2,1),eps));b=-(Js.'*r);
    dz=A\b;dp=scale.*dz;ratio=max(abs(dp)./max(limit,eps));
    if ratio>1,dp=dp/ratio;end
    q=p+dp;value=objective(q);
    if isfinite(value)&&value<old,pNew=q;accepted=true;return;end
    lambda=lambda*10;
end
end

function r=fixed_pair_residual(z,eo,theta,V,x,lib,sensorId,cfg)
u=z(3)*sin(eo(1)*theta)+z(4)*cos(eo(1)*theta)+...
    z(5)*sin(eo(2)*theta)+z(6)*cos(eo(2)*theta);
c=struct('g',z(1),'dx',z(2),'coef',z(3:6));penalty=10*max(std(V),1e-3);
if ~candidate_admissible(c,x,u,cfg),r=penalty*ones(size(V));return;end
r=V-eval_gap_template(lib,z(1),x-z(2)-u,sensorId);r(~isfinite(r))=penalty;
end

function [fit,VFit,rmse]=pack_fit(z,eo,theta,V,x,lib,sensorId,fRot,si)
A=[hypot(z(3),z(4)),hypot(z(5),z(6))];
phi=[atan2(z(4),z(3)),atan2(z(6),z(5))];
u=A(1)*sin(eo(1)*theta+phi(1))+A(2)*sin(eo(2)*theta+phi(2));
VFit=eval_gap_template(lib,z(1),x-z(2)-u,sensorId);r=V-VFit;ok=isfinite(r);
rmse=sqrt(mean(r(ok).^2));fit=struct('g',z(1),'dx',z(2),'A',A,'phi',phi,...
    'eo',eo,'f',eo*fRot,'theta',z,'rmse',rmse,'solve_info',si);
end

function a=merge_struct(a,b)
names=fieldnames(b);for i=1:numel(names),a.(names{i})=b.(names{i});end
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end

function ok=candidate_admissible(c,x,u,cfg)
ok=true;
if isfield(cfg,'supportContract')
    q=x-c.dx-u;s=cfg.supportContract.safe_support_mm;
    gs=cfg.supportContract.static_g_support_mm;
    ok=c.g>=gs(1)&&c.g<=gs(2)&&all(q>=s(1)&q<=s(2));
end
end

function assert_certified_map(x,cfg)
if ~isfield(cfg,'supportContract'),return;end
C=cfg.supportContract;tol=10*eps(max(abs(C.high_certified_mm)));
if any(x<C.high_certified_mm(1)-tol|x>C.high_certified_mm(2)+tol)
    error('support:UncertifiedHighSample','High-speed map contains samples outside its certified domain.');
end
end
