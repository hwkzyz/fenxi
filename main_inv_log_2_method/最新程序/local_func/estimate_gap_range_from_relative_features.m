function gapRange = estimate_gap_range_from_relative_features(highMap, pathModel, cfg)
%ESTIMATE_GAP_RANGE_FROM_RELATIVE_FEATURES  Gap range from high-low features.
% Features are compared as changes relative to the low-speed template sampled
% at the same spatial locations.  This removes much of the window and
% coordinate dependence before using area and width as auxiliary constraints.

if nargin < 3, cfg = struct(); end
xAll = highMap.x_v(:); vAll = highMap.V_a(:);
revAll = highMap.rev_v(:); senAll = highMap.S_v(:);
ok = isfinite(xAll) & isfinite(vAll) & isfinite(revAll) & isfinite(senAll);
xAll = xAll(ok); vAll = vAll(ok); revAll = revAll(ok); senAll = senAll(ok);
[~,~,gid] = unique([revAll,senAll], 'rows', 'stable');
nGroup = max(gid);
xGroup = cell(nGroup,1); vGroup = cell(nGroup,1);
for k=1:nGroup
    ii = gid==k; xGroup{k}=xAll(ii); vGroup{k}=vAll(ii);
end
xTemplate = pathModel.xTemplate(:);
tLow = pathModel.templateLow(:);
g0 = pathModel.pathCal.g0;
dMin = get_cfg(cfg,'gapMin',min(pathModel.templateLib.gapTrain)-g0);
dMax = get_cfg(cfg,'gapMax',max(pathModel.templateLib.gapTrain)-g0);
dGrid = linspace(dMin,dMax,get_cfg(cfg,'gapGridN',41));
shiftGrid = linspace(-get_cfg(cfg,'shiftHalfWidth',0.5), ...
    get_cfg(cfg,'shiftHalfWidth',0.5), get_cfg(cfg,'shiftGridN',11));

% Estimate feature scales from low-speed feature magnitudes.  Height is given
% the largest weight; area and width stabilize the range only when reliable.
baseRef = waveform_features(xTemplate,tLow);
scale = [max(0.03*abs(baseRef(1)),0.005), ...
         max(0.03*abs(baseRef(2)),0.005), ...
         max(0.05*abs(baseRef(3)),0.01)];
weights = get_cfg(cfg,'featureWeights',[1.0,0.35,0.15]);

cost = nan(numel(dGrid),1);
bestShift = nan(numel(dGrid),nGroup);
for ig=1:numel(dGrid)
    tQuery = eval_path_increment_template(pathModel,dGrid(ig),xTemplate);
    groupCost = nan(nGroup,1);
    for k=1:nGroup
        xk=xGroup{k}; vk=vGroup{k};
        if numel(xk)<10, continue; end
        shiftCost=nan(size(shiftGrid));
        for is=1:numel(shiftGrid)
            xq=xk-shiftGrid(is);
            base=interp1(xTemplate,tLow,xq,'pchip',NaN);
            query=interp1(xTemplate,tQuery,xq,'pchip',NaN);
            fObs=waveform_features(xk,vk);
            fBase=waveform_features(xk,base);
            fQuery=waveform_features(xk,query);
            dObs=fObs-fBase; dPred=fQuery-fBase;
            shiftCost(is)=sqrt(sum(weights.*((dObs-dPred)./scale).^2));
        end
        [groupCost(k),bestShift(ig,k)]=min(shiftCost);
        if isfinite(groupCost(k)), bestShift(ig,k)=shiftGrid(bestShift(ig,k)); end
    end
    cost(ig)=median(groupCost,'omitnan');
end
[minCost,iBest]=min(cost);
tol=get_cfg(cfg,'costTolerance',2.5);
inside=cost<=minCost+tol;
dLo=min(dGrid(inside)); dHi=max(dGrid(inside)); dCenter=dGrid(iBest);
half=max([dCenter-dLo,dHi-dCenter,get_cfg(cfg,'minHalfWidth',0.10)]);
gapRange=struct('delta_center',dCenter,'delta_half_width',half, ...
    'delta_lower',max(dMin,dCenter-half),'delta_upper',min(dMax,dCenter+half), ...
    'g_center',g0+dCenter,'g_lower',g0+max(dMin,dCenter-half), ...
    'g_upper',g0+min(dMax,dCenter+half),'delta_grid',dGrid(:), ...
    'cost',cost(:),'best_shift',bestShift,'feature_scale',scale, ...
    'feature_weights',weights,'method',"relative_height_area_width");
end

function f=waveform_features(x,v)
x=x(:);v=v(:);ok=isfinite(x)&isfinite(v);x=x(ok);v=v(ok);
if numel(x)<10,f=[NaN NaN NaN];return;end
[x,ord]=sort(x);v=v(ord);[x,keep]=unique(x,'stable');v=v(keep);
b=prctile(v,5);vp=max(v-b,0);[pk,ip]=max(vp); area=trapz(x,vp);
if pk<=eps,f=[0 0 NaN];return;end
half=pk/2; il=find(vp(1:ip)<=half,1,'last'); ir=find(vp(ip:end)<=half,1,'first');
if isempty(il)||isempty(ir),width=x(end)-x(1);else
    ir=ip+ir-1; width=max(x(ir)-x(il),eps); end
f=[pk area width];
end

function v=get_cfg(s,n,d)
if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
