function low=build_low_speed_templates_binned(dataLow,alpha_k,R_tip,domain,xGrid,statistic,opts)
%BUILD_LOW_SPEED_TEMPLATES_BINNED Controlled spatial-bin template estimator.
% Only the within-bin statistic changes; bin width, coverage, missing-value
% handling and Savitzky-Golay smoothing remain identical.
if nargin<7||isempty(opts),opts=struct();end
if nargin<6||isempty(statistic),statistic='mean';end
dx=get_opt(opts,'gridSpacingMm',.02);minCount=get_opt(opts,'minBinCount',5);span=get_opt(opts,'smoothSpan',9);minCoverage=get_opt(opts,'minCoverage',.90);
if mod(span,2)==0,span=span+1;end
mapRatio=get_opt(opts,'mapWindowRatio',.48);
mapped=map_highspeed_to_space(dataLow,alpha_k,R_tip,domain,0,mapRatio);xOut=xGrid(:);xBin=(domain(1):dx:domain(2)).';if xBin(end)<domain(2)-eps,xBin(end+1)=domain(2);end
sensorIds=unique(mapped.S_v(:),'stable');nS=numel(sensorIds);nB=numel(xBin);Tb=nan(nB,nS);Cb=zeros(nB,nS);Sb=nan(nB,nS);
for is=1:nS
    idx=mapped.S_v==sensorIds(is);x=mapped.x_v(idx);v=mapped.V_a(idx);ib=round((x-xBin(1))/dx)+1;ok=isfinite(x)&isfinite(v)&ib>=1&ib<=nB;ib=ib(ok);v=v(ok);
    for j=1:nB
        q=v(ib==j);Cb(j,is)=numel(q);if numel(q)<minCount,continue;end
        switch lower(char(statistic))
            case 'mean',Tb(j,is)=mean(q);
            case 'median',Tb(j,is)=median(q);
            case 'huber',Tb(j,is)=huber_location(q);
            otherwise,error('Unknown bin statistic: %s',statistic);
        end
        Sb(j,is)=std(q);
    end
    valid=isfinite(Tb(:,is));if nnz(valid)<minCoverage*nB,error('Sensor %d has insufficient spatial-bin coverage.',sensorIds(is));end
    Tb(:,is)=fillmissing(Tb(:,is),'linear','EndValues','nearest');Sb(:,is)=fillmissing(Sb(:,is),'nearest');
    if span>=3,Tb(:,is)=smoothdata(Tb(:,is),'sgolay',span,'omitnan');end
end
T=interp1(xBin,Tb,xOut,'pchip','extrap');C=interp1(xBin,Cb,xOut,'nearest','extrap');S=interp1(xBin,Sb,xOut,'nearest','extrap');
low=struct();low.xGrid=xOut;low.templateBySensor=T;low.templateLow=mean(T,2,'omitnan');low.sensorIds=sensorIds;low.countBySensor=C;low.stdBySensor=S;low.countLow=sum(C,2);low.stdLow=mean(S,2,'omitnan');low.mapped=mapped;low.method=['bin_',char(statistic),'_sg'];low.smooth_span=span;low.grid_spacing_mm=dx;low.map_window_ratio=mapRatio;
end

function mu=huber_location(x)
x=x(isfinite(x));mu=median(x);s=1.4826*median(abs(x-mu));if ~isfinite(s)||s<eps,mu=mean(x);return;end
c=1.345*s;
for it=1:20
    r=x-mu;w=ones(size(r));far=abs(r)>c;w(far)=c./abs(r(far));muNew=sum(w.*x)/sum(w);if abs(muNew-mu)<1e-12,break;end;mu=muNew;
end
end

function v=get_opt(s,n,d)
if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
