function R=run_sync_async_voltage_vp(map,lib,cfg,staticState)
%RUN_SYNC_ASYNC_VOLTAGE_VP One integer-order and one continuous-frequency solver.
t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);fRot=cfg.RPM_high/60;
if isfield(map,'theta_v')&&numel(map.theta_v)==numel(t)
    theta=map.theta_v(:);angleSource="measured_revolution_timing";
else
    theta=2*pi*fRot*t;angleSource="nominal_constant_speed";
end
range=get_field(cfg,'route30DualFrequencyRangeHz',[300 1500]);
eoGrid=ceil(range(1)/fRot):floor(range(2)/fRot);
fixedEO=get_field(cfg,'syncAsyncCandidateEO',[]);
if ~isempty(fixedEO)
    fixedEO=unique(round(fixedEO(:).'));
    if any(~ismember(fixedEO,eoGrid))
        error('syncAsyncCandidateEO must contain integer EO values within the search range.');
    end
    eoGrid=fixedEO;
end
fStep=get_field(cfg,'syncAsyncFrequencyStepHz',10);fGrid=range(1):fStep:range(2);
SS=sin(theta*eoGrid);CS=cos(theta*eoGrid);SA=sin(2*pi*t*fGrid);CA=cos(2*pi*t*fGrid);
gHalf=get_field(cfg,'route30GapHalfWidthMm',.70);gN=get_field(cfg,'syncAsyncGapCount',15);
if isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement
    gGrid=linspace(staticState.gHat-gHalf,staticState.gHat+gHalf,gN);
else
    lo=max(.05,min(lib.gapTrain)-cfg.rawGapSearchMargin);hi=max(lib.gapTrain)+cfg.rawGapSearchMargin;
    gGrid=linspace(max(lo,staticState.gHat-gHalf),min(hi,staticState.gHat+gHalf),gN);
end
bank=repmat(empty_candidate(),0,1);keepPerGap=max(1,round(get_field(cfg,'syncAsyncKeepPerGap',100)));
for g=gGrid
    F0=eval_gap_template(lib,g,x);Fx=eval_gap_derivative(lib,g,x);y=V-F0;
    rows=scan_gap(g,Fx,y,SS,CS,SA,CA,eoGrid,fGrid,fRot,...
        get_field(cfg,'syncAsyncMinSeparationHz',20));
    [~,ord]=sort([rows.sse],'ascend');bank=[bank;rows(ord(1:min(keepPerGap,numel(ord))))]; %#ok<AGROW>
end
[~,ord]=sort([bank.sse],'ascend');iterBudget=max(1,round(get_field(cfg,'syncAsyncIteratedReplayCount',2100)));
iterCount=min(iterBudget,numel(ord));
for ii=1:iterCount
    bank(ord(ii))=iterate_candidate(bank(ord(ii)),theta,t,x,V,lib,...
        get_field(cfg,'syncAsyncIteratedReplayMaxIter',2));
end
[~,ord]=sort([bank.sse],'ascend');refineBudget=max(1,round(get_field(cfg,'syncAsyncRefineCount',50)));
refineCount=min(refineBudget,numel(ord));
best=[];bestRmse=Inf;refinedFrequency=zeros(0,2);refinedRmse=[];refinedFit=cell(0,1);
for ii=1:refineCount
    c=bank(ord(ii));z0=[c.g,c.dx,c.coef(:).',c.f_async];
    gh=get_field(cfg,'route30SingleGapHalfWidthMm',.12);dh=get_field(cfg,'route30DxHalfWidthMm',.5);
    au=get_field(cfg,'route30AmplitudeUpperMm',1.5);fh=max(fStep,get_field(cfg,'syncAsyncRefineHalfWidthHz',15));
    lb=[max(.05,c.g-gh),-dh,-au*ones(1,4),max(range(1),c.f_async-fh)];
    ub=[c.g+gh,dh,au*ones(1,4),min(range(2),c.f_async+fh)];
    if ~(isfield(lib,'fixedPathIncrement')&&lib.fixedPathIncrement)
        lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);
        ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
    end
    [z,si]=solve_lsq_bounded(@(v)mixed_residual(v,c.eo,theta,t,V,x,lib),...
        z0,lb,ub,get_field(cfg,'syncAsyncMaxIter',200),1e-12,1e-12);
    [fit,VFit,rmse]=pack_fit(z,c.eo,theta,t,V,x,lib,fRot,si);
    refinedFrequency(end+1,:)=fit.f;refinedRmse(end+1,1)=rmse; %#ok<AGROW>
    fit.VFit=VFit;refinedFit{end+1,1}=fit; %#ok<AGROW>
    if rmse<bestRmse,best=fit;best.VFit=VFit;bestRmse=rmse;end
end
amplitudeFloor=get_field(cfg,'structuredComponentAmplitudeFloorMm',.075);
eligible=cellfun(@(q)min(q.A)>=amplitudeFloor,refinedFit);
if any(eligible)
    idx=find(eligible);
    [~,j]=min(refinedRmse(idx));best=refinedFit{idx(j)};bestRmse=best.rmse;
end
R=struct('method',"sync_async_voltage_vp",'model_order',2,'mode',"dual_sync_async",...
    'g_used',best.g,'dx_used',best.dx,'A_id',best.A,'phi_id',best.phi,...
    'f_id',best.f,'sync_eo',best.eo,'async_frequency_hz',best.f(2),...
    'rmse',best.rmse,'VFit',best.VFit,'fit',best,'angle_source',angleSource,...
    'candidate_count',numel(bank),'used_sample_count',numel(t),...
    'all_trusted_samples_used',true,'staticState',staticState,...
    'component_amplitude_floor_mm',amplitudeFloor);
C=frequency_solution_confidence(best.f,best.rmse,refinedFrequency(eligible,:),...
    refinedRmse(eligible),numel(t),cfg);
R=merge_struct(R,C);
end

function a=merge_struct(a,b)
names=fieldnames(b);for i=1:numel(names),a.(names{i})=b.(names{i});end
end

function rows=scan_gap(g,Fx,y,SS,CS,SA,CA,eoGrid,fGrid,fRot,minSep)
valid=isfinite(Fx)&isfinite(y);q=-Fx(valid);yv=y(valid);w=q.^2;fy=q.*yv;
SS=SS(valid,:);CS=CS(valid,:);SA=SA(valid,:);CA=CA(valid,:);
g11=sum(w);g1ss=w.'*SS;g1cs=w.'*CS;g1sa=w.'*SA;g1ca=w.'*CA;
ss2=sum(w.*SS.^2,1);sscs=sum(w.*SS.*CS,1);cs2=sum(w.*CS.^2,1);
sa2=sum(w.*SA.^2,1);saca=sum(w.*SA.*CA,1);ca2=sum(w.*CA.^2,1);
sssa=SS.'*(w.*SA);ssca=SS.'*(w.*CA);cssa=CS.'*(w.*SA);csca=CS.'*(w.*CA);
r1=sum(fy);rss=fy.'*SS;rcs=fy.'*CS;rsa=fy.'*SA;rca=fy.'*CA;
rows=repmat(empty_candidate(),0,1);
for ie=1:numel(eoGrid)
    fs=eoGrid(ie)*fRot;
    for jf=1:numel(fGrid)
        if abs(fGrid(jf)-fs)<minSep,continue;end
        G=[g11,g1ss(ie),g1cs(ie),g1sa(jf),g1ca(jf);...
          g1ss(ie),ss2(ie),sscs(ie),sssa(ie,jf),ssca(ie,jf);...
          g1cs(ie),sscs(ie),cs2(ie),cssa(ie,jf),csca(ie,jf);...
          g1sa(jf),sssa(ie,jf),cssa(ie,jf),sa2(jf),saca(jf);...
          g1ca(jf),ssca(ie,jf),csca(ie,jf),saca(jf),ca2(jf)];
        rhs=[r1;rss(ie);rcs(ie);rsa(jf);rca(jf)];
        if any(~isfinite(G),'all')||any(~isfinite(rhs))||rcond(G)<1e-12
            continue;
        end
        beta=G\rhs;
        if any(~isfinite(beta)),continue;end
        sse=max(dot(yv,yv)-dot(rhs,beta),0)/numel(yv);
        rows(end+1,1)=struct('g',g,'dx',beta(1),'eo',eoGrid(ie),...
            'f_sync',fs,'f_async',fGrid(jf),'coef',beta(2:5),'sse',sse); %#ok<AGROW>
    end
end
end

function c=empty_candidate()
c=struct('g',NaN,'dx',NaN,'eo',NaN,'f_sync',NaN,'f_async',NaN,'coef',nan(4,1),'sse',Inf);
end

function c=iterate_candidate(c,theta,t,x,V,lib,maxIter)
dx=c.dx;coef=c.coef(:);S1=sin(c.eo*theta);C1=cos(c.eo*theta);
S2=sin(2*pi*c.f_async*t);C2=cos(2*pi*c.f_async*t);
for it=1:maxIter
    u=coef(1)*S1+coef(2)*C1+coef(3)*S2+coef(4)*C2;
    z=x-dx-u;F0=eval_gap_template(lib,c.g,z);Fx=eval_gap_derivative(lib,c.g,z);
    valid=isfinite(V)&isfinite(F0)&isfinite(Fx);q=-Fx(valid);y=V(valid)-F0(valid);
    X=[q,q.*S1(valid),q.*C1(valid),q.*S2(valid),q.*C2(valid)];
    G=X.'*X;rhs=X.'*y;
    if any(~isfinite(G),'all')||any(~isfinite(rhs))||rcond(G)<1e-12,break;end
    d=G\rhs;
    if any(~isfinite(d)),break;end
    dx=dx+d(1);coef=coef+d(2:5);
end
u=coef(1)*S1+coef(2)*C1+coef(3)*S2+coef(4)*C2;
VFit=eval_gap_template(lib,c.g,x-dx-u);r=V-VFit;good=isfinite(r);
c.dx=dx;c.coef=coef;c.sse=mean(r(good).^2);
end

function r=mixed_residual(z,eo,theta,t,V,x,lib)
u=z(3)*sin(eo*theta)+z(4)*cos(eo*theta)+...
  z(5)*sin(2*pi*z(7)*t)+z(6)*cos(2*pi*z(7)*t);
r=V-eval_gap_template(lib,z(1),x-z(2)-u);r(~isfinite(r))=10*max(std(V),1e-3);
end

function [fit,VFit,rmse]=pack_fit(z,eo,theta,t,V,x,lib,fRot,si)
A=[hypot(z(3),z(4)),hypot(z(5),z(6))];phi=[atan2(z(4),z(3)),atan2(z(6),z(5))];
u=A(1)*sin(eo*theta+phi(1))+A(2)*sin(2*pi*z(7)*t+phi(2));
VFit=eval_gap_template(lib,z(1),x-z(2)-u);r=V-VFit;good=isfinite(r);rmse=sqrt(mean(r(good).^2));
fit=struct('g',z(1),'dx',z(2),'A',A,'phi',phi,'eo',eo,...
    'f',[eo*fRot,z(7)],'theta',z,'rmse',rmse,'solve_info',si);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
