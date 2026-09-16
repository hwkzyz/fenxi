function xBySensor = R5_MakeStaticAuditQueries(models, n)
%R5_MAKESTATICAUDITQUERIES Make a common in-support x grid per model.
if nargin<2 || isempty(n), n=101; end
xBySensor=cell(1,numel(models));
for k=1:numel(models)
    if isfield(models(k),'templateDomainMm') && numel(models(k).templateDomainMm)>=2
        d=double(models(k).templateDomainMm); 
    elseif isfield(models(k),'template') && isfield(models(k).template,'x_grid')
        d=[min(double(models(k).template.x_grid(:))), max(double(models(k).template.x_grid(:)))];
    else
        error('R5:AuditInput','Model %d lacks templateDomainMm.',models(k).sensorId);
    end
    lo=d(1); hi=d(2);
    if isfield(models(k),'responseDomainMm') && isfield(models(k),'coord') && ...
            all(isfinite(models(k).responseDomainMm))
        c=models(k).coord; s=double(c.responseScale); assert(isfinite(s)&&s~=0,'R5:InvalidResponseScale');
        q=models(k).responseDomainMm([1 2])./s + double(c.responseTauMm)+double(c.responseOffsetMm);
        lo=max(lo,min(q)); hi=min(hi,max(q));
    end
    assert(isfinite(lo)&&isfinite(hi)&&hi>lo,'R5:EmptyAuditSupport','No common audit support for sensor %d.',models(k).sensorId);
    xTrial=linspace(lo,hi,max(n,401)).';
    % Leave a small response-coordinate margin for the central derivative audit;
    % the endpoint itself is valid for replay but not for a two-sided probe.
    if isfield(models(k),'coord') && isfield(models(k),'responseDomainMm')
        xMargin=2e-4*max(abs(double(models(k).coord.responseScale)),1);
        if hi-lo>2*xMargin
            xTrial=xTrial(xTrial>lo+xMargin & xTrial<hi-xMargin);
        end
    end
    % overshoot. This avoids declaring a calibration failure merely because
    % the template domain is wider than the response-surface intersection.
    try
        [~,info]=models(k).evaluate(0,xTrial);
        keep=true(size(xTrial));
        if isfield(info,'overshootXmm'), keep=keep & double(info.overshootXmm(:))<=1e-12; end
        if isfield(info,'overshootGapMm'), keep=keep & double(info.overshootGapMm(:))<=1e-12; end
        if isfield(info,'overshootResponseXmm'), keep=keep & double(info.overshootResponseXmm(:))<=1e-12; end
        if isfield(models(k),'isGapSensor') && models(k).isGapSensor
            h=1e-4;
            [yp,~]=models(k).evaluate(h,xTrial);
            [ym,~]=models(k).evaluate(-h,xTrial);
            keep=keep & isfinite(double(yp(:))) & isfinite(double(ym(:)));
        end
        xTrial=xTrial(keep);
    catch
        xTrial=[];
    end
    if numel(xTrial)<max(5,ceil(n/10))
        error('R5:EmptyAuditSupport','No common in-support audit support for sensor %d.',models(k).sensorId);
    end
    pick=round(linspace(1,numel(xTrial),n));
    xBySensor{k}=xTrial(unique(pick));
end
end
