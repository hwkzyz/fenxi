function low = build_low_speed_templates_step04(dataLow,alpha_k,R_tip,domain,xGrid,opts)
%BUILD_LOW_SPEED_TEMPLATES_STEP04 Build independent Step04-like sensor templates.
% Low-speed points are spatially binned, median aggregated, coverage gated,
% and smoothed before they become a displacement-sensitive reference.
if nargin<6,opts=struct();end
dx=get_opt(opts,'gridSpacingMm',median(diff(xGrid(:))));
minCount=get_opt(opts,'minBinCount',5);
smoothSpan=get_opt(opts,'smoothSpan',9);
minCoverage=get_opt(opts,'minCoverage',.90);
mapped=map_highspeed_to_space(dataLow,alpha_k,R_tip,domain,0);
xGrid=xGrid(:);xBuild=(domain(1):dx:domain(2)).';
if xBuild(end)<domain(2)-eps,xBuild(end+1,1)=domain(2);end
sensorIds=unique(mapped.S_v(:),'stable');nBuild=numel(xBuild);nS=numel(sensorIds);
Tbuild=nan(nBuild,nS);Cbuild=zeros(nBuild,nS);Sigbuild=nan(nBuild,nS);kept=zeros(nS,1);
for is=1:nS
    sid=sensorIds(is);idx=mapped.S_v==sid;x=mapped.x_v(idx);v=mapped.V_a(idx);
    ib=round((x-xBuild(1))/dx)+1;
    valid=isfinite(x)&isfinite(v)&ib>=1&ib<=nBuild;ib=ib(valid);v=v(valid);
    for j=1:nBuild
        q=v(ib==j);Cbuild(j,is)=numel(q);
        if numel(q)>=minCount
            Tbuild(j,is)=median(q);
            Sigbuild(j,is)=1.4826*median(abs(q-Tbuild(j,is)));
        end
    end
    valid=isfinite(Tbuild(:,is));
    if nnz(valid)<minCoverage*nBuild
        error('Sensor %d has insufficient Step04-like template coverage.',sid);
    end
    Tbuild(:,is)=fillmissing(Tbuild(:,is),'linear','EndValues','nearest');
    Sigbuild(:,is)=fillmissing(Sigbuild(:,is),'nearest');
    if smoothSpan>1
        span=min(smoothSpan,nBuild);if mod(span,2)==0,span=span-1;end
        Tbuild(:,is)=smoothdata(Tbuild(:,is),'movmedian',span);
        Tbuild(:,is)=smoothdata(Tbuild(:,is),'movmean',span);
    end
    kept(is)=nnz(valid);
end
T=interp1(xBuild,Tbuild,xGrid,'pchip','extrap');
C=interp1(xBuild,Cbuild,xGrid,'nearest','extrap');
Sigma=interp1(xBuild,Sigbuild,xGrid,'nearest','extrap');
low=struct('xGrid',xGrid,'templateBySensor',T,'sensorIds',sensorIds, ...
    'countBySensor',C,'stdBySensor',Sigma,'mapped',mapped, ...
    'num_turns',numel(unique(mapped.rev_v)),'min_bin_count',minCount, ...
    'smooth_span',smoothSpan,'valid_bins_by_sensor',kept);
% Pooled diagnostic only. The forward model uses templateBySensor.
low.templateLow=mean(T,2,'omitnan');low.countLow=sum(C,2);low.stdLow=mean(Sigma,2,'omitnan');
end

function value=get_opt(opts,name,defaultValue)
if isfield(opts,name)&&~isempty(opts.(name)),value=opts.(name);else,value=defaultValue;end
end
