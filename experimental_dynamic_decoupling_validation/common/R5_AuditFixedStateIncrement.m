function audit = R5_AuditFixedStateIncrement(models, xBySensor, dgBySensor, outputFile)
%R5_AUDITFIXEDSTATEINCREMENT Audit one fixed R5 static state without fitting.
% The audit reports template, zero-dg prediction, and finite-dg increment in
% the worker's canonical mV convention.  It is intentionally optimizer-free.
assert(isstruct(models) && ~isempty(models),'R5:MissingModels');
assert(iscell(xBySensor) && numel(xBySensor)==numel(models), ...
    'R5:XContract','xBySensor must match models.');
if nargin<3 || isempty(dgBySensor), dgBySensor=zeros(1,numel(models)); end
dgBySensor=double(dgBySensor(:).'); assert(numel(dgBySensor)==numel(models));
rows=repmat(struct('sensor_id',NaN,'n',0,'template_min_mv',NaN, ...
    'template_max_mv',NaN,'zero_dg_max_abs_mv',NaN,'increment_min_mv',NaN, ...
    'increment_max_mv',NaN,'finite_prediction_min_mv',NaN, ...
    'finite_prediction_max_mv',NaN,'valid_fraction',NaN),numel(models),1);
for k=1:numel(models)
    x=double(xBySensor{k}(:)); assert(~isempty(x),'R5:EmptyAuditGrid');
    yt=models(k).templateEvaluate(x);
    [y0,i0]=models(k).evaluate(0,x); %#ok<ASGLU>
    [y1,i1]=models(k).evaluate(dgBySensor(k),x); %#ok<ASGLU>
    rows(k).sensor_id=double(models(k).sensorId); rows(k).n=numel(x);
    rows(k).template_min_mv=min(yt,[],'omitnan'); rows(k).template_max_mv=max(yt,[],'omitnan');
    rows(k).zero_dg_max_abs_mv=max(abs(y0-yt),[],'omitnan');
    inc=y1-y0; rows(k).increment_min_mv=min(inc,[],'omitnan');
    rows(k).increment_max_mv=max(inc,[],'omitnan');
    rows(k).finite_prediction_min_mv=min(y1,[],'omitnan');
    rows(k).finite_prediction_max_mv=max(y1,[],'omitnan');
    rows(k).valid_fraction=mean(isfinite(y1));
end
audit=struct('schema','R5_FIXED_STATE_INCREMENT_AUDIT_V1', ...
    'rows',rows,'zero_dg_pass',all([rows.zero_dg_max_abs_mv]<1e-8), ...
    'models_sensor_ids',[models.sensorId],'dg_mm',dgBySensor);
if nargin>=4 && ~isempty(outputFile)
    if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
    save(outputFile,'audit','-v7.3');
end
end
