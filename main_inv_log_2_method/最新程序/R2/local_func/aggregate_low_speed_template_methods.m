function low = aggregate_low_speed_template_methods(dataLow,alpha_k,R_tip,domain,xGrid,method,opts)
%AGGREGATE_LOW_SPEED_TEMPLATE_METHODS Compare low-speed template estimators.
% This diagnostic assumes exact spatial/OPR alignment and random noise only.
% Each sensor is kept separate; the method changes only the turn aggregator
% and the final smoothing operation.
if nargin<7 || isempty(opts),opts=struct();end
if nargin<6 || isempty(method),method='mean_raw';end
method=lower(char(method));xGrid=xGrid(:);
smoothSpan=get_opt(opts,'smoothSpan',9);if mod(smoothSpan,2)==0,smoothSpan=smoothSpan+1;end
mapped=map_highspeed_to_space(dataLow,alpha_k,R_tip,domain,0);
sensorIds=unique(mapped.S_v(:),'stable');nS=numel(sensorIds);nX=numel(xGrid);
T=nan(nX,nS);C=zeros(nX,nS);Sg=nan(nX,nS);perTurn=cell(nS,1);
for is=1:nS
    sid=sensorIds(is);rows=unique(mapped.rev_v(mapped.S_v==sid),'stable');Y=nan(numel(rows),nX);
    for ir=1:numel(rows)
        idx=mapped.S_v==sid & mapped.rev_v==rows(ir);
        x=mapped.x_v(idx);v=mapped.V_a(idx);[x,ord]=sort(x(:));v=v(ord);[x,keep]=unique(x,'stable');v=v(keep);
        if numel(x)>=2
            inside=xGrid>=x(1)&xGrid<=x(end);
            Y(ir,inside)=interp1(x,v,xGrid(inside),'pchip');
        end
    end
    perTurn{is}=Y;
    C(:,is)=sum(isfinite(Y),1).';
    switch method
        case {'mean','mean_raw'}
            t=mean(Y,1,'omitnan').';
        case {'median','median_raw','bin_median'}
            t=median(Y,1,'omitnan').';
        case {'mean_sg','mean_sgolay'}
            t=smooth_template(mean(Y,1,'omitnan').',smoothSpan,'sgolay');
        case {'median_sg','median_sgolay'}
            t=smooth_template(median(Y,1,'omitnan').',smoothSpan,'sgolay');
        case {'step04','step04_like'}
            t=smooth_template(smooth_template(median(Y,1,'omitnan').',smoothSpan,'movmedian'),smoothSpan,'movmean');
        otherwise
            error('Unknown low-speed template method: %s',method);
    end
    t=fillmissing(t,'linear','EndValues','nearest');
    T(:,is)=t;
    Sg(:,is)=std(Y,0,1,'omitnan').';
end
low=struct();
low.xGrid=xGrid;low.templateLow=mean(T,2,'omitnan');low.templateBySensor=T;
low.sensorIds=sensorIds;low.countBySensor=C;low.stdBySensor=Sg;low.perTurn=perTurn;
low.mapped=mapped;low.num_turns=numel(unique(mapped.rev_v));low.method=method;low.smooth_span=smoothSpan;
low.countLow=sum(C,2);low.stdLow=mean(Sg,2,'omitnan');
end

function y=smooth_template(y,span,mode)
y=fillmissing(y,'linear','EndValues','nearest');
if span>=3,y=smoothdata(y,mode,span,'omitnan');end
end

function v=get_opt(s,n,d)
if isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
