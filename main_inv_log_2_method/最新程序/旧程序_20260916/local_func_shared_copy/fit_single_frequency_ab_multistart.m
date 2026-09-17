function R=fit_single_frequency_ab_multistart(highMap,lib,cfg,gapState,fSeed,opts)
%FIT_SINGLE_FREQUENCY_AB_MULTISTART  Screen phase/amplitude starts, then fit voltage.
if nargin<6,opts=struct();end
t=highMap.t_v(:);x=highMap.x_v(:);V=highMap.V_a(:);
ampUpper=get_opt(opts,'ampUpperBound',get_cfg(cfg,'finalAmplitudeUpperBoundMm',1.5));
gapHalf=get_opt(opts,'gapHalfWidth',get_cfg(cfg,'finalProjectedJointGapHalfWidth',0.12));
dxHalf=get_opt(opts,'dxHalfWidth',get_cfg(cfg,'finalProjectedJointDxHalfWidth',0.20));
freqHalf=get_opt(opts,'freqHalfWidth',0);
phaseGrid=get_opt(opts,'phaseGrid',(0:23)*pi/12);
ampScales=get_opt(opts,'amplitudeScales',[0.75,0.85,1,1.25,1.4]);
ampSeed=max(get_opt(opts,'amplitudeInit',0.18),get_cfg(cfg,'nonlinearAmpSeedMm',0.18));
physicalAmpMax=min(0.5,ampUpper);
ampStarts=get_opt(opts,'amplitudeStarts',unique([ampSeed*ampScales,0.1:0.1:physicalAmpMax]));
ampStarts=ampStarts(ampStarts>0&ampStarts<=ampUpper);
dxSeeds=get_opt(opts,'dxSeeds',[0,get_opt(opts,'dxInit',0)]);
maxStarts=get_opt(opts,'maxInitialStarts',4);
screenStride=max(1,floor(get_opt(opts,'screenStride',8)));
screenKeep=max(maxStarts,floor(get_opt(opts,'screenKeep',64)));
maxIter=get_opt(opts,'maxIter',300);

if isfield(cfg,'frequencyRangeHz')&&numel(cfg.frequencyRangeHz)>=2
    fDomain=cfg.frequencyRangeHz(1:2);
elseif isfield(cfg,'singleFreqGrid')&&~isempty(cfg.singleFreqGrid)
    fDomain=[min(cfg.singleFreqGrid),max(cfg.singleFreqGrid)];
else
    fDomain=[fSeed,fSeed];
end
lb=[max(0.05,gapState.g_used-gapHalf),-dxHalf,-ampUpper,-ampUpper,max(fDomain(1),fSeed-freqHalf)];
ub=[gapState.g_used+gapHalf,dxHalf,ampUpper,ampUpper,min(fDomain(2),fSeed+freqHalf)];
lb(1)=max(lb(1),min(lib.gapTrain)-cfg.rawGapSearchMargin);
ub(1)=min(ub(1),max(lib.gapTrain)+cfg.rawGapSearchMargin);
penalty=10*max(std(V),1e-3);idx=1:screenStride:numel(t);
res=@(z)single_ab_res(z,t,V,x,lib,penalty);
resScreen=@(z)single_ab_res(z,t(idx),V(idx),x(idx),lib,penalty);

freqStarts=fSeed;
if freqHalf>0
    fineStep=get_opt(opts,'frequencyStartStepHz',0.1);
    fineF=(lb(5):fineStep:ub(5)).';
    if isfield(opts,'inverseMap') && ~isempty(opts.inverseMap)
        inv=opts.inverseMap;
    else
        inv=inverse_map_local_displacement(highMap,lib,gapState.g_used,0,cfg);
    end
    ok=inv.valid(:)&inv.weight(:)>0&isfinite(inv.u_inv(:));
    if nnz(ok)>=6
        tt=t(ok);uu=inv.u_inv(ok);sw=sqrt(inv.weight(ok));fineScore=inf(size(fineF));
        for jf=1:numel(fineF)
            X=[ones(size(tt)),sin(2*pi*fineF(jf)*tt),cos(2*pi*fineF(jf)*tt)];
            b=(X.*sw)\(uu.*sw);r=(uu-X*b).*sw;fineScore(jf)=dot(r,r);
        end
        [~,fo]=sort(fineScore);nFine=get_opt(opts,'frequencyStartTopK',3);
        freqStarts=fineF(fo(1:min(nFine,numel(fo)))).';
    end
end
n=numel(freqStarts)*numel(dxSeeds)*numel(phaseGrid)*numel(ampStarts);
starts=zeros(n,5);score=inf(n,1);k=0;
for jf=1:numel(freqStarts)
for id=1:numel(dxSeeds)
 for ip=1:numel(phaseGrid)
  for ia=1:numel(ampStarts)
   k=k+1;A=ampStarts(ia);ph=phaseGrid(ip);
   z=[gapState.g_used,dxSeeds(id),A*cos(ph),A*sin(ph),freqStarts(jf)];z=min(max(z,lb),ub);
   starts(k,:)=z;r=resScreen(z);score(k)=dot(r,r);
  end
 end
end
end
[~,ord]=sort(score);ord=ord(1:min(screenKeep,numel(ord)));
fullScore=inf(numel(ord),1);
for i=1:numel(ord),r=res(starts(ord(i),:));fullScore(i)=dot(r,r);end
[~,ii]=sort(fullScore);ord=ord(ii(1:min(maxStarts,numel(ii))));
bestSse=inf;theta=[];info=[];
for i=1:numel(ord)
 [z,si]=solve_lsq_bounded(res,starts(ord(i),:),lb,ub,maxIter,1e-12,1e-12);
 r=res(z);sse=dot(r,r);if sse<bestSse,bestSse=sse;theta=z;info=si;end
end
A=hypot(theta(3),theta(4));phi=atan2(theta(4),theta(3));
VFit=eval_gap_template(lib,theta(1),x-theta(2)-A*sin(2*pi*theta(5)*t+phi));
rr=V-VFit;valid=isfinite(rr);
R=struct('g',theta(1),'dx',theta(2),'A',A,'phi',phi,'f',theta(5), ...
 'p',[A,phi,theta(5)],'theta',theta,'VFit',VFit,'rmse',sqrt(mean(rr(valid).^2)), ...
 'invalid_count',sum(~valid),'solve_info',info,'seed_source',"direct_ab_screen");
end

function r=single_ab_res(z,t,V,x,lib,penalty)
A=hypot(z(3),z(4));phi=atan2(z(4),z(3));
y=eval_gap_template(lib,z(1),x-z(2)-A*sin(2*pi*z(5)*t+phi));r=V-y;
bad=~isfinite(r);r(bad)=penalty;
end
function v=get_opt(s,n,d),if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
function v=get_cfg(s,n,d),if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
