function models = R5_Build_SensorConditionedDynamicModels_20241106(sidecarFile, templateFile, targetBlade)
% Build the formal R5 dynamic interface: only dg_s is window-varying.
cfg = Config_20241106();
if nargin < 1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin < 2 || isempty(templateFile), templateFile=cfg.files.lowSpeedTemplate; end
if nargin < 3 || isempty(targetBlade), targetBlade=cfg.case.targetBlade; end
S=load(R5_StageMatForMatlab(sidecarFile,'r5_sidecar'),'SensorConditionedLibrary'); load(R5_StageMatForMatlab(templateFile,'r5_template'),'Template');
L=S.SensorConditionedLibrary; assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
if isfield(Template,'SensorBlade'), templateArray=Template.SensorBlade; else, templateArray=Template.Sensor; end
models=struct('sensorId',{},'bladeId',{},'gapReferenceMm',{},'latent',{}, ...
    'isGapSensor',{},'coord',{},'templateDomainMm',{},'responseDomainMm',{}, ...
    'gapDomainMm',{},'templateBaseline',{},'templateVoltageUnit',{}, ...
    'templateEvaluate',{},'evaluate',{});
analysisSensors=cfg.case.analysisSensors(:).';
gapSensors=cfg.case.gapSensors(:).';
for i=1:numel(L.sensor)
    sc=L.sensor{i}; j=find([sc.state.blade_id]==targetBlade,1);
    if isempty(j), continue; end
    sid=double(sc.sensor_id);
    if isfield(templateArray,'blade_id')
        ib=find([templateArray.sensor_id]==sid & [templateArray.blade_id]==targetBlade,1);
    else
        ib=find([templateArray.sensor_id]==sid,1);
    end
    assert(~isempty(ib),'R5:TemplateMissing','Missing low-speed template for B%d/S%d.',targetBlade,sid);
    tpl=templateArray(ib); reg=sc.registration; st=sc.state(j);
    % Keep the target-blade offset estimated during low-speed localization
    % in the dynamic observation operator.  Omitting it makes the migrated
    % surface and the frozen low-speed template use different x contracts.
    reg.target_x_offset_mm=double(st.target_x_offset_mm);
    coord=struct('templateOffsetMm',0,'responseOffsetMm',reg.target_x_offset_mm, ...
        'responseTauMm',double(reg.tau_mm),'responseScale',double(reg.x_scale), ...
        'responseMuGapPerXMm',get_reg_field(reg,'mu_gap_per_x_mm',0));
    assert(abs(double(tpl.baseline))<1e-12,'R5:TemplateBaselineContract', ...
        'Low-speed templates must be baseline-removed before dynamic fitting.');
    m=struct('sensorId',sid,'bladeId',targetBlade,'gapReferenceMm',st.gap_mm,'isGapSensor',true,...
        'latent',st.z,'coord',coord,'templateDomainMm',[min(tpl.x_grid) max(tpl.x_grid)], ...
        'responseDomainMm',[min(L.x_mm) max(L.x_mm)],'gapDomainMm',[min(L.gap_mm) max(L.gap_mm)], ...
        'templateBaseline',double(tpl.baseline),'templateVoltageUnit',get_template_unit(tpl), ...
        'templateEvaluate',@(xq) eval_template(tpl,xq), ...
        'evaluate',@(dg,xq) eval_one(st.B,L.x_mm,L.gap_mm,L.gref_mm,tpl,reg,st.gap_mm,dg,xq));
    models(end+1)=m; %#ok<AGROW>
end
% Direct/template-only channels remain in the joint waveform objective.  They
% constrain EO, phase, amplitude and dx but do not introduce a spurious gap
% parameter.  Their low-speed template is already the calibrated observation
% model, so the dynamic prediction is simply the template evaluated at xq.
for sid=analysisSensors
    if ismember(sid,gapSensors) || any([models.sensorId]==sid), continue; end
    if isfield(templateArray,'blade_id')
        ib=find([templateArray.sensor_id]==sid & [templateArray.blade_id]==targetBlade,1);
    else
        ib=find([templateArray.sensor_id]==sid,1);
    end
    assert(~isempty(ib),'R5:TemplateMissing','Missing direct low-speed template for B%d/S%d.',targetBlade,sid);
    tpl=templateArray(ib);
    coord=struct('templateOffsetMm',0,'responseOffsetMm',0,'responseTauMm',0,'responseScale',1);
    models(end+1)=struct('sensorId',sid,'bladeId',targetBlade,'gapReferenceMm',NaN,...
        'latent',NaN,'isGapSensor',false,'coord',coord, ...
        'templateDomainMm',[min(tpl.x_grid) max(tpl.x_grid)], ...
        'responseDomainMm',[min(tpl.x_grid) max(tpl.x_grid)],'gapDomainMm',[NaN NaN], ...
        'templateBaseline',double(tpl.baseline),'templateVoltageUnit',get_template_unit(tpl), ...
        'templateEvaluate',@(xq) eval_template(tpl,xq), ...
        'evaluate',@(dg,xq) eval_direct_template(tpl,xq)); %#ok<AGROW>
end
models=sort_models_by_sensor(models,analysisSensors);
assert(~isempty(models),'R5:NoDynamicSensors','No target-blade sensor states found.');
end

function u=get_template_unit(tpl)
u='V'; if isfield(tpl,'voltageUnit') && ~isempty(tpl.voltageUnit), u=char(tpl.voltageUnit); end
end

function y=eval_template(tpl,xq)
y=interp1(double(tpl.x_grid(:)),template_to_mv(tpl),double(xq(:)),'pchip',NaN);
end

function models=sort_models_by_sensor(models,ids)
[~,ord]=ismember(ids,[models.sensorId]);
assert(all(ord>0),'R5:SensorRoleMismatch');
models=models(ord);
end

function [y,info]=eval_direct_template(tpl,xq)
xq=xq(:); xg=double(tpl.x_grid(:));
lo=min(xg); hi=max(xg);
y=interp1(xg,template_to_mv(tpl),xq,'pchip',NaN);
info=struct('xQuery',xq,'xRegistered',xq,'gapLowMm',NaN(size(xq)),...
    'noGapMv',y,'gapIncrementMv',zeros(size(xq)),...
    'gapHighMm',NaN(size(xq)),'gapClippedMm',NaN(size(xq)),...
    'overshootGapMm',zeros(size(xq)),'overshootResponseXmm',zeros(size(xq)), ...
    'overshootXmm',max(lo-xq,0)+max(xq-hi,0));
end

function y=template_to_mv(tpl)
% V1 templates may be stored either in volts or millivolts.  Convert once
% at this boundary so an mV template is never multiplied by 1000 again.
u='V';
if isfield(tpl,'voltageUnit') && ~isempty(tpl.voltageUnit)
    u=lower(strtrim(char(tpl.voltageUnit)));
end
y=double(tpl.v_grid(:));
if any(strcmp(u,{'v','volt','volts'}))
    y=1000*y;
elseif ~any(strcmp(u,{'mv','millivolt','millivolts'}))
    error('R5:UnknownTemplateVoltageUnit','Unsupported template voltage unit: %s',u);
end
end

function [y,info]=eval_one(B,xGrid,gapGrid,gref,tpl,reg,gL,dg,xq)
xq=xq(:); dg=double(dg); assert(isscalar(dg),'R5:ScalarDg');
if isfield(reg,'target_x_offset_mm') && isfinite(reg.target_x_offset_mm)
    xOffset=double(reg.target_x_offset_mm);
else
    xOffset=0;
end
xResponse=xq-xOffset;
xF=double(reg.x_scale).*(xResponse-double(reg.tau_mm));
xLo=min(double(tpl.x_grid)); xHi=max(double(tpl.x_grid)); xLowRaw=xq;
% The measured low-speed template remains in its own OPR coordinate.
% target_offset belongs only to the calibrated response-surface map.
yL=interp1(double(tpl.x_grid(:)),double(tpl.v_grid(:))*1000,xLowRaw,'pchip',NaN);
% The localization contract supplies one sensor-conditioned effective
% reference gap for each blade.  Do not re-introduce the calibration
% tilt path here: localization was performed against the scalar-gap
% surface, so applying mu(x-tau) only in the dynamic stage creates a
% static-template/model mismatch and contaminates dg.  The required
% clearance freedom remains the independent dynamic increment dg.
    mu=get_reg_field(reg,'mu_gap_per_x_mm',0);
    gHL=gL + mu.*(xq-double(reg.tau_mm));
gHH=gHL+dg;
gMin=min(gapGrid); gMax=max(gapGrid);
f0=eval_surface(B,xGrid,gref,gHL,xF); f1=eval_surface(B,xGrid,gref,gHH,xF);
y=yL+double(reg.voltage_gain).*(f1-f0);
valid=xLowRaw>=xLo & xLowRaw<=xHi & xF>=min(xGrid) & xF<=max(xGrid) & ...
    gHL>=gMin & gHL<=gMax & gHH>=gMin & gHH<=gMax;
y(~valid)=NaN;
info=struct('xQuery',xq,'xRegistered',xF,'gapLowMm',gHL,'gapHighMm',gHH,...
    'noGapMv',yL,'gapIncrementMv',double(reg.voltage_gain).*(f1-f0),...
    'gapClippedMm',gHH,'overshootGapMm',max(gMin-gHH,0)+max(gHH-gMax,0),...
    'overshootXmm',max(xLo-xLowRaw,0)+max(xLowRaw-xHi,0), ...
    'overshootResponseXmm',max(min(xGrid)-xF,0)+max(xF-max(xGrid),0));
end

function v=get_reg_field(s,name,d)
v=d; if isfield(s,name) && isfinite(s.(name)), v=double(s.(name)); end
end


function y=eval_surface(B,xGrid,gref,g,xq)
v=zeros(size(xq));
for k=1:3
    bk=interp1(xGrid(:),B(:,k),xq,'linear',NaN);
    if k==1, v=v+bk; elseif k==2, v=v+bk./g; else, v=v+bk.*log(g/gref); end
end
y=v;
end
