function audit = R5_Audit_StaticObservationContract(models, xBySensor, toleranceMv, minFullSupport, allowGapExtrapolation)
%R5_AUDIT_STATICOBSERVATIONCONTRACT Replay the frozen observation operator.
% This is a pre-fit gate: zero gap increment must reproduce the frozen
% low-speed template and every declared query must remain inside support.
if nargin < 3 || isempty(toleranceMv), toleranceMv = 1e-8; end
if nargin < 4 || isempty(minFullSupport), minFullSupport = 0.80; end
if nargin < 5 || isempty(allowGapExtrapolation), allowGapExtrapolation = false; end
assert(isstruct(models) && ~isempty(models),'R5:NoModels');
assert(iscell(xBySensor) && numel(xBySensor)==numel(models), ...
    'R5:AuditInput','xBySensor must contain one query vector per model.');
rows = repmat(struct('sensor_id',NaN,'zero_gap_max_abs_mv',NaN, ...
    'zero_gap_pass',false,'support_pass',false,'metadata_pass',false,'finite_count',0, ...
    'template_baseline',NaN,'template_unit_ok',false,'coordinate_pass',false, ...
    'template_replay_max_abs_mv',NaN,'template_replay_pass',false, ...
    'coordinate_max_abs_mm',NaN,'derivative_pass',false, ...
    'full_template_support_fraction',NaN,'full_template_evaluable_fraction',NaN, ...
    'full_template_support_domain',[NaN NaN]),numel(models),1);
for k=1:numel(models)
    x = double(xBySensor{k}(:));
    [y,info] = models(k).evaluate(0,x);
    y=double(y(:));
    if isfield(info,'noGapMv')
        noGap=double(info.noGapMv(:));
    else
        % Sensor-conditioned models return the frozen template replay as the
        % zero-increment prediction; retain a non-fatal compatibility path
        % while still checking finiteness, support and derivatives.
        noGap=y;
    end
    d = y-noGap; good=isfinite(d);
    rows(k).sensor_id=double(models(k).sensorId);
    hasCommonMetadata=isfield(models(k),'templateDomainMm') && ...
        numel(models(k).templateDomainMm)>=2 && ...
        all(isfinite(double(models(k).templateDomainMm(1:2)))) && ...
        isfield(models(k),'coord') && ...
        all(isfinite([models(k).coord.responseScale, ...
        models(k).coord.responseTauMm,models(k).coord.responseOffsetMm]));
    isGap=isfield(models(k),'isGapSensor') && logical(models(k).isGapSensor);
    hasGapMetadata=~isGap || (isfield(models(k),'gapReferenceMm') && ...
        isfield(models(k),'gapDomainMm') && numel(models(k).gapDomainMm)>=2 && ...
        isfield(models(k),'responseDomainMm') && numel(models(k).responseDomainMm)>=2 && ...
        all(isfinite([double(models(k).gapReferenceMm), ...
        double(models(k).gapDomainMm(1:2)),double(models(k).responseDomainMm(1:2))])) && ...
        double(models(k).gapReferenceMm)>=min(double(models(k).gapDomainMm(1:2))) && ...
        double(models(k).gapReferenceMm)<=max(double(models(k).gapDomainMm(1:2))));
    rows(k).metadata_pass=hasCommonMetadata && hasGapMetadata;
    % Audit the complete frozen template domain separately.  The query grid
    % may intentionally be trimmed to an in-support subset, but that must
    % remain visible: a short overlap is a calibration limitation, not a
    % successful full-domain replay.
    if isfield(models(k),'templateDomainMm') && numel(models(k).templateDomainMm)>=2
        xd=linspace(double(models(k).templateDomainMm(1)),double(models(k).templateDomainMm(2)),401).';
        try
            [yd,id]=models(k).evaluate(0,xd);
            yd=double(yd(:));
            okd=true(size(xd));
            if isfield(id,'overshootXmm'), okd=okd & double(id.overshootXmm(:))<=1e-12; end
            gapSupported=okd;
            if isfield(id,'overshootGapMm'), gapSupported=gapSupported & double(id.overshootGapMm(:))<=1e-12; end
            rows(k).full_template_evaluable_fraction=nnz(isfinite(yd))/numel(yd);
            calibratedDomain=gapSupported & isfinite(yd);
            rows(k).full_template_support_fraction=nnz(calibratedDomain)/numel(calibratedDomain);
            if allowGapExtrapolation
                okd=isfinite(yd);
            else
                okd=calibratedDomain;
            end
            if any(okd), rows(k).full_template_support_domain=[min(xd(okd)) max(xd(okd))]; end
        catch
            rows(k).full_template_support_fraction=0;
        end
    end
    % Unit/baseline metadata are part of the static contract.  A model may
    % still be numerically self-consistent after a hidden 1000x conversion,
    % so this check is deliberately independent of the zero-gap identity.
    if isfield(models(k),'templateBaseline')
        rows(k).template_baseline=double(models(k).templateBaseline);
    elseif isfield(models(k),'template') && isfield(models(k).template,'baseline')
        rows(k).template_baseline=double(models(k).template.baseline);
    end
    rows(k).template_unit_ok=~isfield(models(k),'templateVoltageUnit') || ...
        strcmpi(string(models(k).templateVoltageUnit),'V');
    rows(k).coordinate_pass=isfield(models(k),'coord') && ...
        all(isfinite([models(k).coord.responseScale,models(k).coord.responseTauMm,models(k).coord.responseOffsetMm]));
    if rows(k).coordinate_pass && isfield(info,'xRegistered')
        c=models(k).coord;
        xExpected=double(c.responseScale).*(x-double(c.responseTauMm)-double(c.responseOffsetMm));
        xc=double(info.xRegistered(:));
        if numel(xc)==numel(x) && all(isfinite(xc))
            rows(k).coordinate_max_abs_mm=max(abs(xc-xExpected));
            rows(k).coordinate_pass=rows(k).coordinate_max_abs_mm<=1e-10;
        else
            rows(k).coordinate_pass=false;
        end
    end
    if isfield(models(k),'templateEvaluate') && ~isempty(models(k).templateEvaluate)
        yTemplate=double(models(k).templateEvaluate(x));
        dTemplate=noGap-yTemplate;
        qTemplate=isfinite(dTemplate);
        if any(qTemplate), rows(k).template_replay_max_abs_mv=max(abs(dTemplate(qTemplate))); end
        rows(k).template_replay_pass=all(qTemplate) && ...
            rows(k).template_replay_max_abs_mv<=toleranceMv;
    end
    h=1e-4;
    derivativeMode='central';
    if isfield(models(k),'isGapSensor') && models(k).isGapSensor && ...
            isfield(models(k),'gapReferenceMm') && isfield(models(k),'gapDomainMm') && ...
            all(isfinite([models(k).gapReferenceMm models(k).gapDomainMm]))
        gRef=double(models(k).gapReferenceMm);
        gDom=double(models(k).gapDomainMm(:));
        lowerMargin=gRef-min(gDom); upperMargin=max(gDom)-gRef;
        if lowerMargin < h && upperMargin >= h
            derivativeMode='forward';
        elseif upperMargin < h && lowerMargin >= h
            derivativeMode='backward';
        elseif lowerMargin < h && upperMargin < h
            h=max(min(lowerMargin,upperMargin),eps);
            if upperMargin >= lowerMargin
                derivativeMode='forward';
            else
                derivativeMode='backward';
            end
        end
    end
    try
        y0=double(y(:));
        if strcmp(derivativeMode,'central')
            [yp,ip]=models(k).evaluate(h,x);
            [ym,im]=models(k).evaluate(-h,x);
            yp=double(yp(:)); ym=double(ym(:));
            plusOk=all(isfinite(yp)) && audit_info_inside(ip);
            minusOk=all(isfinite(ym)) && audit_info_inside(im);
            if plusOk && minusOk
                deriv=(yp-ym)/(2*h);
            elseif plusOk
                deriv=(yp-y0)/h;
            elseif minusOk
                deriv=(y0-ym)/h;
            else
                deriv=NaN(size(y0));
            end
        elseif strcmp(derivativeMode,'forward')
            [yp,ip]=models(k).evaluate(h,x); yp=double(yp(:));
            if all(isfinite(yp)) && audit_info_inside(ip), deriv=(yp-y0)/h; else, deriv=NaN(size(y0)); end
        else
            [ym,im]=models(k).evaluate(-h,x); ym=double(ym(:));
            if all(isfinite(ym)) && audit_info_inside(im), deriv=(y0-ym)/h; else, deriv=NaN(size(y0)); end
        end
        rows(k).derivative_pass=all(isfinite(deriv));
    catch
        rows(k).derivative_pass=false;
    end
    rows(k).finite_count=nnz(good);
    if any(good), rows(k).zero_gap_max_abs_mv=max(abs(d(good))); end
    rows(k).zero_gap_pass=isfinite(rows(k).zero_gap_max_abs_mv) && ...
        rows(k).zero_gap_max_abs_mv <= toleranceMv;
    if isfield(info,'overshootXmm') && isfield(info,'overshootGapMm')
        if allowGapExtrapolation
            % x is intentionally clamped to the finite response-coordinate
            % interval, while positive-gap extrapolation is evaluated by the
            % analytic basis.  The contract is therefore finite/evaluable;
            % both overshoot vectors remain available for reporting.
            support=all(isfinite(y));
        else
            support=all(info.overshootXmm==0) && all(info.overshootGapMm==0);
        end
        if isfield(info,'overshootResponseXmm'), support=support && all(info.overshootResponseXmm==0); end
        rows(k).support_pass=support;
    else
        rows(k).support_pass=all(isfinite(y));
    end
end
audit=struct('schema','R5_STATIC_OBSERVATION_CONTRACT_AUDIT_V2', ...
    'toleranceMv',toleranceMv,'rows',rows, ...
    'minimumFullTemplateSupportFraction',minFullSupport, ...
    'allowGapExtrapolation',logical(allowGapExtrapolation), ...
    'pass',all([rows.zero_gap_pass] & [rows.support_pass] & ...
        [rows.metadata_pass] & [rows.template_unit_ok] & [rows.coordinate_pass] & ...
        [rows.template_replay_pass] & [rows.derivative_pass] & ...
        ((allowGapExtrapolation .* [rows.full_template_evaluable_fraction] + ...
          (~allowGapExtrapolation) .* [rows.full_template_support_fraction]) >= minFullSupport)));
end

function tf=audit_info_inside(info)
tf=true;
for name={'overshootXmm','overshootGapMm','overshootResponseXmm'}
    if isfield(info,name{1})
        v=double(info.(name{1})(:));
        tf=tf && all(isfinite(v)) && all(v<=1e-12);
    end
end
end
