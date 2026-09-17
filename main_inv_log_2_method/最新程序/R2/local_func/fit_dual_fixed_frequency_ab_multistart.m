function R=fit_dual_fixed_frequency_ab_multistart(highMap,lib,cfg,gapState,fPair,seed,opts)
%FIT_DUAL_FIXED_FREQUENCY_AB_MULTISTART  Fixed-frequency dual voltage fit.
if nargin<7,opts=struct();end
t=highMap.t_v(:);x=highMap.x_v(:);V=highMap.V_a(:);f=sort(fPair(:).');
ampUpper=get_opt(opts,'ampUpperBound',get_cfg(cfg,'finalAmplitudeUpperBoundMm',1.5));
gapHalf=get_opt(opts,'gapHalfWidth',get_cfg(cfg,'finalProjectedJointGapHalfWidth',0.12));
dxHalf=get_opt(opts,'dxHalfWidth',get_cfg(cfg,'finalProjectedJointDxHalfWidth',0.20));
phaseStarts=get_opt(opts,'phaseStarts',(0:23)*pi/12);
dxSeeds=0;
maxStarts=get_opt(opts,'maxInitialStarts',8);screenStride=get_opt(opts,'screenStride',16);
screenKeep=max(maxStarts,get_opt(opts,'screenKeep',96));maxIter=get_opt(opts,'maxIter',400);
lb=[max(0.05,gapState.g_used-gapHalf),-dxHalf,-ampUpper,-ampUpper,-ampUpper,-ampUpper];
ub=[gapState.g_used+gapHalf,dxHalf,ampUpper,ampUpper,ampUpper,ampUpper];
lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
penalty=10*max(std(V),1e-3);res=@(z)dual_res(z,t,V,x,lib,f,penalty);ii=1:screenStride:numel(t);
resScreen=@(z)dual_res(z,t(ii),V(ii),x(ii),lib,f,penalty);
A0=max(seed.A(:).',[0.10,0.10]);physical=0.1:0.05:min(0.5,ampUpper);
Astart1=unique([A0(1),physical]);Astart2=unique([A0(2),physical]);
n=numel(dxSeeds)*numel(phaseStarts)^2*numel(Astart1)*numel(Astart2);starts=zeros(n,6);score=inf(n,1);k=0;
for id=1:numel(dxSeeds)
 for ip1=1:numel(phaseStarts)
  for ip2=1:numel(phaseStarts)
   for ia1=1:numel(Astart1)
    for ia2=1:numel(Astart2)
     k=k+1;A1=Astart1(ia1);A2=Astart2(ia2);p1=phaseStarts(ip1);p2=phaseStarts(ip2);
     z=[gapState.g_used,dxSeeds(id),A1*cos(p1),A1*sin(p1),A2*cos(p2),A2*sin(p2)];z=min(max(z,lb),ub);
     starts(k,:)=z;r=resScreen(z);score(k)=dot(r,r);
    end
   end
  end
 end
end
[~,ord]=sort(score);ord=ord(1:min(screenKeep,numel(ord)));full=inf(numel(ord),1);
for i=1:numel(ord),r=res(starts(ord(i),:));full(i)=dot(r,r);end
[~,q]=sort(full);ord=ord(q(1:min(maxStarts,numel(q))));best=inf;theta=[];info=[];
for i=1:numel(ord)
 [z,si]=solve_lsq_bounded(res,starts(ord(i),:),lb,ub,maxIter,1e-12,1e-12);r=res(z);s=dot(r,r);
 if s<best,best=s;theta=z;info=si;end
end
A=[hypot(theta(3),theta(4)),hypot(theta(5),theta(6))];ph=[atan2(theta(4),theta(3)),atan2(theta(6),theta(5))];
u=A(1)*sin(2*pi*f(1)*t+ph(1))+A(2)*sin(2*pi*f(2)*t+ph(2));VFit=eval_gap_template(lib,theta(1),x-theta(2)-u);rr=V-VFit;ok=isfinite(rr);
R=struct('g',theta(1),'dx',theta(2),'A',A,'phi',ph,'f',f,'p',[A(1),ph(1),f(1),A(2),ph(2),f(2)], ...
 'theta',theta,'VFit',VFit,'rmse',sqrt(mean(rr(ok).^2)),'invalid_count',sum(~ok),'solve_info',info,'seed_source',"inverse_direct_ab");
end
function r=dual_res(z,t,V,x,lib,f,penalty)
A1=hypot(z(3),z(4));p1=atan2(z(4),z(3));A2=hypot(z(5),z(6));p2=atan2(z(6),z(5));
u=A1*sin(2*pi*f(1)*t+p1)+A2*sin(2*pi*f(2)*t+p2);r=V-eval_gap_template(lib,z(1),x-z(2)-u);bad=~isfinite(r);r(bad)=penalty;
end
function v=get_opt(s,n,d),if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
function v=get_cfg(s,n,d),if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
