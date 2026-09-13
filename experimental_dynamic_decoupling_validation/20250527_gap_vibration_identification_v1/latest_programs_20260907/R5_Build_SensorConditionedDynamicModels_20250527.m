function models = R5_Build_SensorConditionedDynamicModels_20250527(sidecarFile, templateFile, targetBlade, v1Calibration)
% Build the formal R5 dynamic interface: only dg_s is window-varying.
cfg = Config_20250527();
if nargin < 4, v1Calibration = []; end
if nargin < 1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin < 2 || isempty(templateFile), templateFile=cfg.files.lowSpeedTemplate; end
if nargin < 3 || isempty(targetBlade), targetBlade=cfg.case.targetBlade; end
S=load(R5_StageMatForMatlab(sidecarFile,'r5_sidecar'),'SensorConditionedLibrary'); load(R5_StageMatForMatlab(templateFile,'r5_template'),'Template');
L=S.SensorConditionedLibrary; assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
if isfield(Template,'SensorBlade'), templateArray=Template.SensorBlade; else, templateArray=Template.Sensor; end
% All three analysis channels in the 20250527 condition are clearance
% channels.  The shared worker uses this explicit role flag to allocate the
% sensor-wise dynamic dg parameters.
models=struct('sensorId',{},'bladeId',{},'gapReferenceMm',{},'latent',{}, ...
    'isGapSensor',{},'coord',{},'templateDomainMm',{},'responseDomainMm',{}, ...
    'gapDomainMm',{},'templateBaseline',{},'templateVoltageUnit',{}, ...
    'templateEvaluate',{},'evaluate',{});
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
    if ~isempty(v1Calibration) && isfield(v1Calibration,'sensor')
        iv=find([v1Calibration.sensor.sensorId]==sid,1);
        if ~isempty(iv)
            cv=v1Calibration.sensor(iv);
            if isfield(cv,'tauMm'), reg.tau_mm=cv.tauMm; end
            if isfield(cv,'xScale'), reg.x_scale=cv.xScale; end
            if isfield(cv,'muGapPerXMm'), reg.mu_gap_per_x_mm=cv.muGapPerXMm; end
            if isfield(cv,'voltageGain'), reg.voltage_gain=cv.voltageGain; end
        end
    end
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
        'evaluate',@(dg,xq) eval_one(st.B,L.x_mm,L.gap_mm,scalar_value(L.gref_mm),tpl,reg,scalar_value(st.gap_mm),dg,xq));
    models(end+1)=m; %#ok<AGROW>
end

function v=scalar_value(v)
v=double(v(:)); v=v(isfinite(v)); assert(~isempty(v),'R5:MissingScalar'); v=mean(v);
end
assert(~isempty(models),'R5:NoDynamicSensors','No target-blade sensor states found.');
end

function u=get_template_unit(tpl)
u='V'; if isfield(tpl,'voltageUnit') && ~isempty(tpl.voltageUnit), u=char(tpl.voltageUnit); end
end

function y=eval_template(tpl,xq)
y=interp1(double(tpl.x_grid(:)),template_values_mv(tpl),double(xq(:)),'pchip',NaN);
end

function y=template_values_mv(tpl)
u='v'; if isfield(tpl,'voltageUnit') && ~isempty(tpl.voltageUnit), q=tpl.voltageUnit; if iscell(q), q=q{1}; end; u=lower(strtrim(char(q))); end
y=double(tpl.v_grid(:));
if any(strcmp(u,{'v','volt','volts'})), y=1000*y;
elseif ~any(strcmp(u,{'mv','millivolt','millivolts'})), error('R5:UnknownTemplateVoltageUnit','Unsupported template voltage unit: %s',u); end
if isfield(tpl,'baseline'), y=y-double(tpl.baseline); end
end

function [y,info]=eval_one(B,xGrid,gapGrid,gref,tpl,reg,gL,dg,xq)
xq=xq(:); dg=double(dg); assert(isscalar(dg),'R5:ScalarDg');
% The localization stage estimates a fixed target-blade/template offset in
% the sensor-conditioned library coordinate.  It is a frozen static
% registration term, not a dynamic degree of freedom; omitting it here
% makes the migrated surface query a different x-coordinate from the one
% used during low-speed localization.
if isfield(reg,'target_x_offset_mm') && isfinite(reg.target_x_offset_mm)
    xOffset=double(reg.target_x_offset_mm);
else
    xOffset=0;
end
xResponse=xq-xOffset;
xF=double(reg.x_scale(1)).*(xResponse-double(reg.tau_mm(1)));
xLo=min(double(tpl.x_grid)); xHi=max(double(tpl.x_grid)); xLowRaw=xq;
% Localization maps the measured template abscissa xa into the response
% surface with x_scale*(xa-tau-target_offset). The measured template stays
% on xa itself. Applying target_offset to this template shifts it twice.
yL=interp1(double(tpl.x_grid(:)),template_values_mv(tpl),xLowRaw,'pchip',NaN);
% The target low-speed template for this condition retains a measurable
% local clearance slope.  Use that calibrated, frozen baseline when the
% condition contract enables it; only dg is window-varying.
    mu=get_reg_field(reg,'mu_gap_per_x_mm',0); mu=double(mu(1));
gHL=gL + mu.*(xq-double(reg.tau_mm(1)));
gHH=gHL+dg;
gMin=min(gapGrid); gMax=max(gapGrid);
f0=eval_surface(B,xGrid,gref,gHL,xF); f1=eval_surface(B,xGrid,gref,gHH,xF);
y=yL+double(reg.voltage_gain(1)).*(f1-f0);
valid=xLowRaw>=xLo & xLowRaw<=xHi & xF>=min(xGrid) & xF<=max(xGrid) & ...
    gHL>=gMin & gHL<=gMax & gHH>=gMin & gHH<=gMax;
y(~valid)=NaN;
info=struct('xQuery',xq,'xRegistered',xF,'gapLowMm',gHL,'gapHighMm',gHH,...
    'noGapMv',yL,'gapIncrementMv',double(reg.voltage_gain(1)).*(f1-f0),...
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
