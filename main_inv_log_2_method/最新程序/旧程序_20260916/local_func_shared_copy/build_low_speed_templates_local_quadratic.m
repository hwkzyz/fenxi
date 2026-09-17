function low=build_low_speed_templates_local_quadratic(dataLow,alpha_k,R_tip,domain,xGrid,opts)
%BUILD_LOW_SPEED_TEMPLATES_LOCAL_QUADRATIC Support-aware local quadratic template.
% A spatial bin is used only to obtain stable local means. The template and
% derivative are then the intercept and slope of a weighted local quadratic,
% which retains sub-bin spatial phase information.
if nargin<6||isempty(opts),opts=struct();end
dx=get_opt(opts,'gridSpacingMm',.02);
windowMm=get_opt(opts,'windowMm',.12);
minCount=get_opt(opts,'minBinCount',5);
mapRatio=get_opt(opts,'mapWindowRatio',.52);
mapped=map_highspeed_to_space(dataLow,alpha_k,R_tip,domain,0,mapRatio);
xGrid=xGrid(:);xBin=(domain(1):dx:domain(2)).';
if xBin(end)<domain(2)-eps,xBin(end+1)=domain(2);end
sensorIds=unique(mapped.S_v(:),'stable');nS=numel(sensorIds);nB=numel(xBin);
Tb=nan(nB,nS);Cb=zeros(nB,nS);
for is=1:nS
    q=mapped.S_v==sensorIds(is);x=mapped.x_v(q);v=mapped.V_a(q);
    ib=round((x-xBin(1))/dx)+1;ok=isfinite(x)&isfinite(v)&ib>=1&ib<=nB;
    ib=ib(ok);v=v(ok);
    for j=1:nB
        z=v(ib==j);Cb(j,is)=numel(z);
        if numel(z)>=minCount,Tb(j,is)=mean(z);end
    end
    if nnz(isfinite(Tb(:,is)))<0.90*nB
        error('Sensor %d has insufficient spatial-bin coverage.',sensorIds(is));
    end
    Tb(:,is)=fillmissing(Tb(:,is),'linear','EndValues','nearest');
end
T=nan(numel(xGrid),nS);D=T;
for is=1:nS
    for ix=1:numel(xGrid)
        h=xBin-xGrid(ix);use=abs(h)<=windowMm/2;
        if nnz(use)<3
            [~,ord]=sort(abs(h));use=false(size(h));use(ord(1:min(3,numel(h))))=true;
        end
        hs=h(use);ys=Tb(use,is);scale=max(windowMm/2,dx);
        u=hs/scale;w=exp(-0.5*(u/.75).^2);
        X=[ones(size(u)),hs,hs.^2];X=X.*sqrt(w);
        beta=X\(ys.*sqrt(w));T(ix,is)=beta(1);D(ix,is)=beta(2);
    end
end
C=interp1(xBin,Cb,xGrid,'nearest','extrap');
low=struct('xGrid',xGrid,'templateBySensor',T,'derivativeBySensor',D,...
    'templateLow',mean(T,2,'omitnan'),'derivativeLow',mean(D,2,'omitnan'),...
    'sensorIds',sensorIds,'countBySensor',C,'stdBySensor',nan(size(T)),...
    'countLow',sum(C,2),'stdLow',nan(size(xGrid)),'mapped',mapped,...
    'method','local_quadratic','grid_spacing_mm',dx,...
    'local_window_mm',windowMm,'map_window_ratio',mapRatio);
end

function v=get_opt(s,n,d)
if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
