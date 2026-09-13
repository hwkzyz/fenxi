function [p,c] = R5_Evaluate_CanonicalForward(z,f,B,M)
%R5_EVALUATECANONICALFORWARD Unified low-template differential-gap operator.
% z = [A, phase, dx, dg_1 ... dg_S].  Static registration, template,
% baseline-gap path and gain live in M; only z is dynamic.  This function
% is the shared forward backend for experimental R5 and synthetic checks.
eo=double(B.EO); if ~isfinite(eo), eo=f/B.rotFreqMeanHz; end
phase=eo*B.Theta+z(2)+2*pi*(f-eo*B.rotFreqMeanHz)*(B.T-mean(B.T));
u=z(1)*sin(phase); p=nan(size(B.V));
c=struct('foundation',nan(size(B.V)),'noGapCurrent',nan(size(B.V)), ...
    'noGapBase',nan(size(B.V)),'gapIncrement',nan(size(B.V)), ...
    'xCurrent',B.X-z(3)-u,'xBase',nan(size(B.V)),'uCurrent',u);
for i=1:numel(M)
    q=B.sensorIndex==i;
    [p(q),d]=M(i).evaluate(z(3+i),c.xCurrent(q));
    c.noGapCurrent(q)=d.noGapMv; c.gapIncrement(q)=d.gapIncrementMv;
    if isfield(B,'useV1DynamicIncrement') && B.useV1DynamicIncrement
        [c.gapIncrement(q),okv]=v1_dynamic_increment(B,M(i),z(3+i),c.xCurrent(q));
        p(q)=d.noGapMv+c.gapIncrement(q); p(q(~okv))=NaN;
    end
end

function [d,ok]=v1_dynamic_increment(B,m,dg,x)
d=nan(size(x)); ok=false(size(x));
if ~isfield(B,'v1ResponseSurface')||isempty(B.v1ResponseSurface)||~isfield(B,'v1SensorCalibration')||isempty(B.v1SensorCalibration), return; end
R=B.v1ResponseSurface; cal=B.v1SensorCalibration; ii=find([cal.sensorId]==m.sensorId,1); if isempty(ii), return; end
c0=cal(ii); k=c0.xScale; tau=c0.tauMm; mu=c0.muGapPerXMm; gain=c0.voltageGain;
jj=find(B.sensorIds==m.sensorId,1); if isempty(jj), return; end
g0=B.v1G0BySensor(jj); xr=k.*(x-tau); ge=g0+mu.*(x-tau); g1=ge+dg;
d=gain.*(v1_surface(R,g1,xr)-v1_surface(R,ge,xr));
ok=isfinite(d)&xr>=min(R.xGrid)&xr<=max(R.xGrid)&g1>=min(R.gTrainMm)&g1<=max(R.gTrainMm);
end
if isfield(B,'useV1LowTemplate') && B.useV1LowTemplate
    for i=1:numel(M)
        q=B.sensorIndex==i;
        p(q)=v1_low_curve(B,M(i),c.xCurrent(q))+c.gapIncrement(q);
    end
end
if isfield(B,'useNestedFoundationAnchor') && B.useNestedFoundationAnchor
    assert(all(isfinite([B.anchorEO B.anchorF B.anchorA B.anchorPhi B.anchorDx])), ...
        'R5:IncompleteFoundationAnchor');
    z0=[B.anchorA B.anchorPhi B.anchorDx zeros(1,numel(M))];
    eo0=B.anchorEO; f0=B.anchorF;
    u0=z0(1)*sin(eo0*B.Theta+z0(2)+2*pi*(f0-eo0*B.rotFreqMeanHz)*(B.T-mean(B.T)));
    p0=nan(size(B.V));
    for i=1:numel(M)
        q=B.sensorIndex==i;
        [p0(q),~]=M(i).evaluate(0,B.X(q)-z0(3)-u0(q));
    end
    q=isfinite(B.anchorV)&isfinite(p)&isfinite(p0);
    p(q)=B.anchorV(q)+p(q)-p0(q);
    c.foundation=B.anchorV; c.noGapBase=p0;
else
    z0=z; z0(4:end)=0; u0=z0(1)*sin(phase-z(2));
    for i=1:numel(M)
        q=B.sensorIndex==i;
        [c.noGapBase(q),~]=M(i).evaluate(0,B.X(q)-z0(3)-u0(q));
    end
end

function v=v1_low_curve(B,m,x)
T=B.v1Template; SI=B.v1SensorInfo;
v=nan(size(x)); x=x(:);
if ~isempty(SI) && isstruct(SI) && isfield(SI,'sensorId')
    ii=find([SI.sensorId]==m.sensorId,1);
    if ~isempty(ii) && isfield(SI(ii),'directBaseX') && isfield(SI(ii),'directBaseF0')
        xx=double(SI(ii).directBaseX(:)); yy=double(SI(ii).directBaseF0(:));
        if numel(xx)>=8 && numel(xx)==numel(yy)
        v=interp1(xx,yy,min(max(x,min(xx)),max(xx)),'pchip',nan);
        end
    end
end
if isstruct(T) && (isfield(T,'SensorBlade') || isfield(T,'Sensor'))
    if isfield(T,'SensorBlade'), arr=T.SensorBlade; else, arr=T.Sensor; end
    ii=find([arr.sensor_id]==m.sensorId,1);
    if ~isempty(ii)
        a=arr(ii); xx=double(a.x_grid(:));
        v=interp1(xx,template_values_mv(a),min(max(x,min(xx)),max(xx)),'pchip',nan);
    end
end

function y=template_values_mv(a)
u='v';
if isfield(a,'voltageUnit') && ~isempty(a.voltageUnit)
    q=a.voltageUnit; if iscell(q), q=q{1}; end
    u=lower(strtrim(char(q)));
end
y=double(a.v_grid(:));
if any(strcmp(u,{'v','volt','volts'})), y=1000*y;
elseif ~any(strcmp(u,{'mv','millivolt','millivolts'})), error('R5:UnknownTemplateVoltageUnit','Unsupported template voltage unit: %s',u); end
if isfield(a,'baseline'), y=y-double(a.baseline); end
end
% A baseline-subtracted V1 template is the no-gap curve.  The old
% gap_only contract also contained a frozen operating-gap increment
% F(g0)-F(gref).  Add that term here, while keeping dg exclusively in the
% R5 increment above; otherwise migration changes the absolute voltage
% baseline by hundreds of mV and corrupts f/EO/A.
if ~isempty(B.v1ResponseSurface) && ~isempty(B.v1SensorCalibration)
    cal=B.v1SensorCalibration; ii=[];
    if isstruct(cal) && isfield(cal,'sensorId')
        ii=find([cal.sensorId]==m.sensorId,1);
    end
    if ~isempty(ii)
        c0=cal(ii); R=B.v1ResponseSurface;
        if isfield(c0,'xScale'), k=double(c0.xScale); else, k=1; end
        if isfield(c0,'tauMm'), tau=double(c0.tauMm); else, tau=0; end
        if isfield(c0,'muGapPerXMm'), mu=double(c0.muGapPerXMm); else, mu=0; end
        if isfield(c0,'voltageGain'), gain=double(c0.voltageGain); else, gain=1; end
        if ~isempty(B.v1G0BySensor)
            jj=find(B.sensorIds==m.sensorId,1);
            if ~isempty(jj) && isfinite(B.v1G0BySensor(jj)), g0=B.v1G0BySensor(jj); else, g0=double(c0.g0Mm); end
        else
            g0=double(c0.g0Mm);
        end
        xr=k.*(x-tau); ge=g0+mu.*(x-tau);
        gref=double(R.g0Mm);
        f0=v1_surface(R,ge,xr); fr=v1_surface(R,gref+zeros(size(xr)),xr);
        v=v+gain.*(f0-fr);
    end
end
end

function y=v1_surface(R,g,x)
% Legacy response surfaces use B0+B1/g+B2*log(g/gref) in mV.
X=double(R.xGrid(:)); C=double(R.coeff);
b0=interp1(X,C(:,1),x,'linear',NaN);
b1=interp1(X,C(:,2),x,'linear',NaN);
b2=interp1(X,C(:,3),x,'linear',NaN);
y=b0+b1./g+b2.*log(g./double(R.g0Mm));
end
end
