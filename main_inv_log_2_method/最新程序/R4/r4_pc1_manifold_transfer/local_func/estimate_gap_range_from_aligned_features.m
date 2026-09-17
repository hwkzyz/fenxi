function gapRange = estimate_gap_range_from_aligned_features(highMap, pathModel, cfg)
%ESTIMATE_GAP_RANGE_FROM_ALIGNED_FEATURES  Align first, then use multi-features.

if nargin < 3, cfg = struct(); end
xAll=highMap.x_v(:); vAll=highMap.V_a(:);
revAll=highMap.rev_v(:); senAll=highMap.S_v(:);
ok=isfinite(xAll)&isfinite(vAll)&isfinite(revAll)&isfinite(senAll);
xAll=xAll(ok);vAll=vAll(ok);revAll=revAll(ok);senAll=senAll(ok);
[~,~,gid]=unique([revAll,senAll],'rows','stable'); nGroup=max(gid);
xTemplate=pathModel.xTemplate(:); tLow=pathModel.templateLow(:);
shiftHalf=get_cfg(cfg,'shiftHalfWidth',0.5); shiftGrid=linspace(-shiftHalf,shiftHalf,get_cfg(cfg,'shiftGridN',41));
xGroup=cell(nGroup,1); vGroup=cell(nGroup,1); shiftHat=nan(nGroup,1);
obsDelta=nan(nGroup,3);
for k=1:nGroup
    ii=gid==k; xk=xAll(ii);vk=vAll(ii); xGroup{k}=xk;vGroup{k}=vk;
    if numel(xk)<10,continue;end
    score=nan(size(shiftGrid));
    for is=1:numel(shiftGrid)
        base=interp1(xTemplate,tLow,xk-shiftGrid(is),'pchip',NaN);
        good=isfinite(base)&isfinite(vk);
        if nnz(good)>8
            a=vk(good)-mean(vk(good)); b=base(good)-mean(base(good));
            score(is)=(a'*b)/(norm(a)*norm(b)+eps);
        end
    end
    [~,ib]=max(score); shiftHat(k)=shiftGrid(ib);
    base=interp1(xTemplate,tLow,xk-shiftHat(k),'pchip',NaN);
    obsDelta(k,:)=waveform_features(xk,vk)-waveform_features(xk,base);
end
keep=isfinite(shiftHat)&all(isfinite(obsDelta),2); xGroup=xGroup(keep);vGroup=vGroup(keep);shiftHat=shiftHat(keep);obsDelta=obsDelta(keep,:);
if size(obsDelta,1)<2,error('Too few valid groups after coarse alignment.');end

baseRef=waveform_features(xTemplate,tLow);
scale=[max(0.03*abs(baseRef(1)),0.005),max(0.03*abs(baseRef(2)),0.005),max(0.05*abs(baseRef(3)),0.01)];
weights=get_cfg(cfg,'featureWeights',[1,0.35,0.15]);
g0=pathModel.pathCal.g0; dMin=get_cfg(cfg,'gapMin',min(pathModel.templateLib.gapTrain)-g0); dMax=get_cfg(cfg,'gapMax',max(pathModel.templateLib.gapTrain)-g0);
dGrid=linspace(dMin,dMax,get_cfg(cfg,'gapGridN',41)); cost=nan(size(dGrid));
for ig=1:numel(dGrid)
    tq=eval_path_increment_template(pathModel,dGrid(ig),xTemplate); gc=nan(size(shiftHat));
    for k=1:numel(shiftHat)
        xk=xGroup{k}; base=interp1(xTemplate,tLow,xk-shiftHat(k),'pchip',NaN); query=interp1(xTemplate,tq,xk-shiftHat(k),'pchip',NaN);
        predDelta=waveform_features(xk,query)-waveform_features(xk,base);
        gc(k)=sqrt(sum(weights.*((obsDelta(k,:)-predDelta)./scale).^2));
    end
    cost(ig)=median(gc,'omitnan');
end
[minCost,ib]=min(cost); inside=cost<=minCost+get_cfg(cfg,'costTolerance',2.5); dLo=min(dGrid(inside));dHi=max(dGrid(inside));dc=dGrid(ib); half=max([dc-dLo,dHi-dc,get_cfg(cfg,'minHalfWidth',0.10)]);
gapRange=struct('delta_center',dc,'delta_half_width',half,'delta_lower',max(dMin,dc-half),'delta_upper',min(dMax,dc+half),'g_center',g0+dc,'g_lower',g0+max(dMin,dc-half),'g_upper',g0+min(dMax,dc+half),'delta_grid',dGrid(:),'cost',cost(:),'shift_hat',shiftHat,'obs_delta_features',obsDelta,'feature_scale',scale,'feature_weights',weights,'method',"aligned_relative_height_area_width");
end

function f=waveform_features(x,v)
x=x(:);v=v(:);ok=isfinite(x)&isfinite(v);x=x(ok);v=v(ok);if numel(x)<10,f=[NaN NaN NaN];return;end
[x,o]=sort(x);v=v(o);[x,keep]=unique(x,'stable');v=v(keep);b=prctile(v,5);vp=max(v-b,0);[pk,ip]=max(vp);ar=trapz(x,vp);if pk<=eps,f=[0 0 NaN];return;end
hh=pk/2;il=find(vp(1:ip)<=hh,1,'last');ir=find(vp(ip:end)<=hh,1,'first');if isempty(il)||isempty(ir),wd=x(end)-x(1);else,ir=ip+ir-1;wd=max(x(ir)-x(il),eps);end;f=[pk ar wd];
end

function v=get_cfg(s,n,d),if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end,end
