function R=run_single_sync_voltage_vp(map,lib,cfg,staticState)
%RUN_SINGLE_SYNC_VOLTAGE_VP Integer-order single-frequency voltage solver.
% All trusted samples are retained. The synchronous basis uses measured
% shaft angle when available and falls back to nominal constant-speed angle.
t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);
fRot=cfg.RPM_high/60;
if isfield(map,'theta_v')&&numel(map.theta_v)==numel(t)
    theta=map.theta_v(:);
    angleSource="measured_revolution_timing";
else
    theta=2*pi*fRot*t;
    angleSource="nominal_constant_speed";
end
range=get_field(cfg,'route30SingleFrequencyRangeHz',[5 1500]);
eoGrid=ceil(range(1)/fRot):floor(range(2)/fRot);
if isempty(eoGrid),error('invlog2:EmptySynchronousOrderGrid','No integer EO is inside the requested frequency range.');end
fixedEO=get_field(cfg,'singleSyncCandidateEO',[]);
if ~isempty(fixedEO)
    fixedEO=unique(round(fixedEO(:).'));
    if any(~ismember(fixedEO,eoGrid))
        error('singleSyncCandidateEO must contain integer EO values within the search range.');
    end
    eoGrid=fixedEO;
end
S=sin(theta*eoGrid);C=cos(theta*eoGrid);
gCenter=staticState.gHat;gapHalf=get_field(cfg,'route30GapHalfWidthMm',.70);
if isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement
    gapGrid=linspace(gCenter-gapHalf,gCenter+gapHalf,get_field(cfg,'route30GapCount',61));
else
    lo=max(.05,min(lib.gapTrain)-cfg.rawGapSearchMargin);
    hi=max(lib.gapTrain)+cfg.rawGapSearchMargin;
    gapGrid=linspace(max(lo,gCenter-gapHalf),min(hi,gCenter+gapHalf),...
        get_field(cfg,'route30GapCount',61));
end
bank=repmat(struct('g',NaN,'eo',NaN,'f',NaN,'dx',NaN,'a',NaN,'b',NaN,'sse',Inf),0,1);
keepPerGap=get_field(cfg,'syncSingleKeepPerGap',numel(eoGrid));
for ig=1:numel(gapGrid)
    g=gapGrid(ig);F0=eval_gap_template(lib,g,x);Fx=eval_gap_derivative(lib,g,x);y=V-F0;
    valid=isfinite(Fx)&isfinite(y)&abs(Fx)>eps;
    if nnz(valid)<8,continue;end
    % Use the same ordinary voltage-domain least-squares metric as the
    % complete model. Sensitivity weighting is useful for diagnostics, but
    % Fx^4 weighting can rank large-amplitude nonlinear excursions poorly.
    q=-Fx(valid);yv=y(valid);d=q.^2;
    Sv=S(valid,:);Cv=C(valid,:);nf=numel(eoGrid);G=zeros(3,3,nf);
    G(1,1,:)=sum(d);G(1,2,:)=d.'*Sv;G(2,1,:)=G(1,2,:);
    G(1,3,:)=d.'*Cv;G(3,1,:)=G(1,3,:);G(2,2,:)=sum(d.*Sv.^2,1);
    G(2,3,:)=sum(d.*Sv.*Cv,1);G(3,2,:)=G(2,3,:);G(3,3,:)=sum(d.*Cv.^2,1);
    rs=q.*yv;rhs=[repmat(sum(rs),1,nf);rs.'*Sv;rs.'*Cv];
    beta=reshape(pagemldivide(G,reshape(rhs,3,1,nf)),3,nf);
    sse=max(sum(yv.^2)-sum(rhs.*beta,1),0);sse(~isfinite(sse))=Inf;
    [~,ord]=sort(sse,'ascend');ord=ord(1:min(keepPerGap,numel(ord)));
    for j=ord
        bank(end+1,1)=struct('g',g,'eo',eoGrid(j),'f',eoGrid(j)*fRot,...
            'dx',beta(1,j),'a',beta(2,j),'b',beta(3,j),'sse',sse(j)); %#ok<AGROW>
    end
end
if isempty(bank),error('invlog2:NoSynchronousCandidate','No valid synchronous candidate was found.');end
[bank,directReplayCount]=direct_sync_replay(bank,eoGrid,theta,x,V,lib,cfg);
[~,ord]=sort([bank.sse],'ascend');refineCount=min(get_field(cfg,'syncSingleRefineCount',3),numel(ord));
iterCount=min(get_field(cfg,'syncSingleIteratedReplayCount',numel(ord)),numel(ord));
for ii=1:iterCount
    bank(ord(ii))=iterated_sync_candidate(bank(ord(ii)),theta,x,V,lib,...
        get_field(cfg,'syncSingleIteratedReplayMaxIter',2));
end
[~,ord]=sort([bank.sse],'ascend');
best=[];bestRmse=Inf;refinedFrequency=[];refinedRmse=[];
for k=1:refineCount
    c=bank(ord(k));z0=[c.g,c.dx,c.a,c.b];
    gh=get_field(cfg,'route30SingleGapHalfWidthMm',.12);dh=get_field(cfg,'route30DxHalfWidthMm',.5);
    au=get_field(cfg,'route30AmplitudeUpperMm',1.5);
    lb=[max(.05,c.g-gh),-dh,-au,-au];ub=[c.g+gh,dh,au,au];
    if ~(isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement)
        lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);
        ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
    end
    residual=@(z) sync_residual(z,c.eo,theta,V,x,lib);
    [z,si]=solve_lsq_bounded(residual,z0,lb,ub,get_field(cfg,'route30SingleMaxIter',160),1e-12,1e-12);
    A=hypot(z(3),z(4));phi=atan2(z(4),z(3));u=A*sin(c.eo*theta+phi);
    VFit=eval_gap_template(lib,z(1),x-z(2)-u);r=V-VFit;good=isfinite(r);
    rmse=sqrt(mean(r(good).^2));
    refinedFrequency(end+1,1)=c.f;refinedRmse(end+1,1)=rmse; %#ok<AGROW>
    if rmse<bestRmse
        bestRmse=rmse;best=struct('g',z(1),'dx',z(2),'A',A,'phi',phi,...
            'eo',c.eo,'f',c.f,'theta',z,'VFit',VFit,'rmse',rmse,'solve_info',si);
    end
end

R=struct('method',"single_sync_voltage_vp",'model_order',1,'mode',"single_sync",...
    'g_used',best.g,'dx_used',best.dx,'A_id',best.A,'phi_id',best.phi,...
    'f_id',best.f,'eo_id',best.eo,'rmse',best.rmse,'VFit',best.VFit,...
    'fit',best,'angle_source',angleSource,'eo_grid',eoGrid(:),...
    'candidate_count',numel(bank),'direct_replay_count',directReplayCount,...
    'used_sample_count',numel(t),...
    'all_trusted_samples_used',true,'staticState',staticState);
C=frequency_solution_confidence(best.f,best.rmse,refinedFrequency,refinedRmse,numel(t),cfg);
R=merge_struct(R,C);
end

function a=merge_struct(a,b)
names=fieldnames(b);for i=1:numel(names),a.(names{i})=b.(names{i});end
end

function [out,nReplay]=direct_sync_replay(bank,eoGrid,theta,x,V,lib,cfg)
% Every integer order receives exact-forward amplitude/phase replay at
% several independently ranked gap candidates. This avoids using a small-
% displacement VP score as the final gate for a large-amplitude response.
amplitudeGrid=get_field(cfg,'syncSingleAmplitudeSeedGridMm',[.1 .2 .3 .4 .5]);
phaseGrid=get_field(cfg,'syncSinglePhaseSeedGrid',0:pi/4:(2*pi-pi/4));
out=bank;nReplay=0;[AA,PP]=ndgrid(amplitudeGrid,phaseGrid);
avec=AA(:).';pvec=PP(:).';
for j=1:numel(bank)
    c=bank(j);U=sin(c.eo*theta+pvec).*avec;
    VFit=eval_gap_template(lib,c.g,x-c.dx-U);r=V-VFit;
    valid=isfinite(r);r(~valid)=0;sse=sum(r.^2,1)./max(sum(valid,1),1);
    [bestSse,k]=min(sse);bestA=avec(k);bestPhi=pvec(k);nReplay=nReplay+numel(sse);
    c.a=bestA*cos(bestPhi);c.b=bestA*sin(bestPhi);c.sse=bestSse;out(j)=c;
end
end

function c=iterated_sync_candidate(c,theta,x,V,lib,maxIter)
dx=c.dx;coef=[c.a;c.b];S=sin(c.eo*theta);C=cos(c.eo*theta);
for it=1:maxIter
    u=coef(1)*S+coef(2)*C;z=x-dx-u;
    F0=eval_gap_template(lib,c.g,z);Fx=eval_gap_derivative(lib,c.g,z);
    valid=isfinite(V)&isfinite(F0)&isfinite(Fx);q=-Fx(valid);y=V(valid)-F0(valid);
    X=[q,q.*S(valid),q.*C(valid)];delta=X\y;
    dx=dx+delta(1);coef=coef+delta(2:3);
end
VFit=eval_gap_template(lib,c.g,x-dx-coef(1)*S-coef(2)*C);
r=V-VFit;good=isfinite(r);c.sse=mean(r(good).^2);c.dx=dx;c.a=coef(1);c.b=coef(2);
end

function r=sync_residual(z,eo,theta,V,x,lib)
A=hypot(z(3),z(4));phi=atan2(z(4),z(3));
r=V-eval_gap_template(lib,z(1),x-z(2)-A*sin(eo*theta+phi));
r(~isfinite(r))=10*max(std(V),1e-3);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
