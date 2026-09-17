function models = R5_Build_SensorConditionedDynamicModels_20250527(sidecarFile, templateFile, targetBlade)
% Build the formal R5 dynamic interface: only dg_s is window-varying.
cfg = Config_20250527();
if nargin < 1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin < 2 || isempty(templateFile), templateFile=cfg.files.lowSpeedTemplate; end
if nargin < 3 || isempty(targetBlade), targetBlade=cfg.case.targetBlade; end
S=load(R5_StageMatForMatlab(sidecarFile,'r5_sidecar'),'SensorConditionedLibrary'); load(R5_StageMatForMatlab(templateFile,'r5_template'),'Template');
L=S.SensorConditionedLibrary; assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
assert(isfield(L,'staticTiltCarrier') && isfield(L,'staticTiltFrozen') && ...
    isfield(L,'spaceGapSlopeUsed') && isfield(L,'allowDynamicTilt'), ...
    'R5:StaticTiltContract','Sidecar lacks the static tilt contract.');
assert(~logical(L.allowDynamicTilt), 'R5:DynamicTiltForbidden', ...
    'Dynamic tilt is forbidden in the unified adapter.');
assert(strcmpi(char(L.staticTiltCarrier),'absolute_surface_state') && ...
    ~logical(L.spaceGapSlopeUsed), 'R5:StaticTiltCarrier', ...
    'Dynamic adapter requires an absolute frozen surface; explicit gap slope is forbidden.');
% All three analysis channels in the 20250527 condition are clearance
% channels.  The shared worker uses this explicit role flag to allocate the
% sensor-wise dynamic dg parameters.
models=struct('sensorId',{},'bladeId',{},'gapReferenceMm',{},'latent',{},'isGapSensor',{}, ...
    'templateDomainMm',{},'responseDomainMm',{},'gapDomainMm',{},'templateVoltageUnit',{}, ...
    'coord',{},'templateEvaluate',{},'evaluate',{});
for i=1:numel(L.sensor)
    sc=L.sensor{i}; j=find([sc.state.blade_id]==targetBlade,1);
    if isempty(j), continue; end
    sid=double(sc.sensor_id); ib=find([Template.SensorBlade.sensor_id]==sid & ...
        [Template.SensorBlade.blade_id]==targetBlade,1);
    assert(~isempty(ib),'R5:TemplateMissing','Missing low-speed template for B%d/S%d.',targetBlade,sid);
    tpl=Template.SensorBlade(ib); reg=sc.registration; st=sc.state(j);
    if isfield(st,'target_x_offset_mm') && cfg.r5.useTargetXOffsetInDynamic
        reg.target_x_offset_mm=double(st.target_x_offset_mm);
    else
        reg.target_x_offset_mm=0;
    end
    m=struct('sensorId',sid,'bladeId',targetBlade,'gapReferenceMm',st.gap_mm,'isGapSensor',true,...
        'templateDomainMm',[min(double(tpl.x_grid(:))) max(double(tpl.x_grid(:)))], ...
        'responseDomainMm',[min(double(L.x_mm(:))) max(double(L.x_mm(:)))], ...
        'gapDomainMm',[min(double(L.gap_mm(:))) max(double(L.gap_mm(:)))], ...
        'templateVoltageUnit','V','coord',struct('responseScale',double(reg.x_scale), ...
        'responseTauMm',double(reg.tau_mm),'responseOffsetMm',double(reg.target_x_offset_mm)), ...
        'templateEvaluate',@(xq) interp1(double(tpl.x_grid(:)),double(tpl.v_grid(:))*1000, ...
        min(max(double(xq(:))-double(reg.target_x_offset_mm),min(double(tpl.x_grid(:)))),max(double(tpl.x_grid(:)))), ...
        'pchip',NaN), ...
        'latent',st.z,'evaluate',@(dg,xq) eval_one(st.B,L.x_mm,L.gap_mm,L.gref_mm,...
        tpl,reg,st.gap_mm,dg,xq));
    models(end+1)=m; %#ok<AGROW>
end
assert(~isempty(models),'R5:NoDynamicSensors','No target-blade sensor states found.');
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
xF=double(reg.x_scale).*(xResponse-double(reg.tau_mm));
xLo=min(double(tpl.x_grid)); xHi=max(double(tpl.x_grid)); xLowRaw=xq-xOffset; xClip=min(max(xLowRaw,xLo),xHi);
% The frozen target offset maps the dynamic OPR coordinate into the target
% blade's low-speed/template coordinate and the response-surface coordinate;
% both terms must use this same registration.
% CoreBundlePreview.v is the absolute sensor voltage.  Keep the template in
% that same voltage frame.  Subtracting tpl.baseline here alone introduced
% a sensor-wise constant error (59.8--89.1 mV in window 1), which the joint
% fit then incorrectly absorbed into A, dx and dg_s.
yL=interp1(double(tpl.x_grid(:)),double(tpl.v_grid(:))*1000,xClip,'pchip',NaN);
% The target low-speed template for this condition retains a measurable
% local clearance slope.  Use that calibrated, frozen baseline when the
% condition contract enables it; only dg is window-varying.
% Static spatial variation is carried by the localized absolute surface;
% dynamic evaluation starts from scalar g0 and never adds explicit mu/tau.
gHL=repmat(gL,size(xq));
gHH=gHL+dg;
gMin=min(gapGrid); gMax=max(gapGrid);
% Keep the complete waveform in the dynamic objective.  The surface is only
% measured on [gMin,gMax], so use the V1 boundary extension for the fit, but
% retain an explicit mask and magnitude for every clipped sample.  This keeps
% the objective complete without presenting the extension as new calibration
% data.
g0Eval=min(max(gHL,gMin),gMax); g1Eval=min(max(gHH,gMin),gMax);
f0=eval_surface(B,xGrid,gref,g0Eval,xF); f1=eval_surface(B,xGrid,gref,g1Eval,xF);
y=yL+double(reg.voltage_gain).*(f1-f0);
noGapMv=yL;
gapOvershoot=max(max(gMin-gHL,0)+max(gHL-gMax,0), ...
    max(gMin-gHH,0)+max(gHH-gMax,0));
gapIncrementMv=y-noGapMv;
info=struct('xQuery',xq,'xRegistered',xF,'gapLowMm',gHL,'gapHighMm',gHH,...
    'gapEvaluatedMm',g1Eval,'overshootGapMm',gapOvershoot,...
    'extrapolatedGapMask',(gHL<gMin | gHL>gMax | gHH<gMin | gHH>gMax),...
    'overshootXmm',max(xLo-xLowRaw,0)+max(xLowRaw-xHi,0),...
    'noGapMv',noGapMv,'gapIncrementMv',gapIncrementMv);
end


function y=eval_surface(B,xGrid,gref,g,xq)
v=zeros(size(xq));
for k=1:3
    bk=interp1(xGrid(:),B(:,k),min(max(xq,min(xGrid)),max(xGrid)),'linear',NaN);
    if k==1, v=v+bk; elseif k==2, v=v+bk./g; else, v=v+bk.*log(g/gref); end
end
y=v;
end
