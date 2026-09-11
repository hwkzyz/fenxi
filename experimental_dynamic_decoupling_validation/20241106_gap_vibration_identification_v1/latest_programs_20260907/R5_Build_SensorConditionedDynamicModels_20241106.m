function models = R5_Build_SensorConditionedDynamicModels_20241106(sidecarFile, templateFile, targetBlade)
% Build the formal R5 dynamic interface: only dg_s is window-varying.
cfg = Config_20241106();
if nargin < 1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin < 2 || isempty(templateFile), templateFile=cfg.files.lowSpeedTemplate; end
if nargin < 3 || isempty(targetBlade), targetBlade=cfg.case.targetBlade; end
S=load(R5_StageMatForMatlab(sidecarFile,'r5_sidecar'),'SensorConditionedLibrary'); load(R5_StageMatForMatlab(templateFile,'r5_template'),'Template');
L=S.SensorConditionedLibrary; assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
models=struct('sensorId',{},'bladeId',{},'gapReferenceMm',{},'latent',{},'isGapSensor',{},'evaluate',{});
analysisSensors=cfg.case.analysisSensors(:).';
gapSensors=cfg.case.gapSensors(:).';
for i=1:numel(L.sensor)
    sc=L.sensor{i}; j=find([sc.state.blade_id]==targetBlade,1);
    if isempty(j), continue; end
    sid=double(sc.sensor_id); ib=find([Template.SensorBlade.sensor_id]==sid & ...
        [Template.SensorBlade.blade_id]==targetBlade,1);
    assert(~isempty(ib),'R5:TemplateMissing','Missing low-speed template for B%d/S%d.',targetBlade,sid);
    tpl=Template.SensorBlade(ib); reg=sc.registration; st=sc.state(j);
    % Keep the target-blade offset estimated during low-speed localization
    % in the dynamic observation operator.  Omitting it makes the migrated
    % surface and the frozen low-speed template use different x contracts.
    if isfield(st,'target_x_offset_mm') && cfg.r5.useTargetXOffsetInDynamic
        reg.target_x_offset_mm=double(st.target_x_offset_mm);
    else
        reg.target_x_offset_mm=0;
    end
    m=struct('sensorId',sid,'bladeId',targetBlade,'gapReferenceMm',st.gap_mm,'isGapSensor',true,...
        'latent',st.z,'evaluate',@(dg,xq) eval_one(st.B,L.x_mm,L.gap_mm,L.gref_mm,...
        tpl,reg,st.gap_mm,dg,xq,cfg.r5.useStaticGapSlopeInDynamic,st.mu_gap_per_x_mm,st.gap_tau_mm));
    models(end+1)=m; %#ok<AGROW>
end
% Direct/template-only channels remain in the joint waveform objective.  They
% constrain EO, phase, amplitude and dx but do not introduce a spurious gap
% parameter.  Their low-speed template is already the calibrated observation
% model, so the dynamic prediction is simply the template evaluated at xq.
for sid=analysisSensors
    if ismember(sid,gapSensors) || any([models.sensorId]==sid), continue; end
    ib=find([Template.SensorBlade.sensor_id]==sid & ...
        [Template.SensorBlade.blade_id]==targetBlade,1);
    assert(~isempty(ib),'R5:TemplateMissing','Missing direct low-speed template for B%d/S%d.',targetBlade,sid);
    tpl=Template.SensorBlade(ib);
    models(end+1)=struct('sensorId',sid,'bladeId',targetBlade,'gapReferenceMm',NaN,...
        'latent',NaN,'isGapSensor',false,'evaluate',@(dg,xq) eval_direct_template(tpl,xq)); %#ok<AGROW>
end
models=sort_models_by_sensor(models,analysisSensors);
assert(~isempty(models),'R5:NoDynamicSensors','No target-blade sensor states found.');
end

function models=sort_models_by_sensor(models,ids)
[~,ord]=ismember(ids,[models.sensorId]);
assert(all(ord>0),'R5:SensorRoleMismatch');
models=models(ord);
end

function [y,info]=eval_direct_template(tpl,xq)
xq=xq(:); xg=double(tpl.x_grid(:));
lo=min(xg); hi=max(xg); xc=min(max(xq,lo),hi);
y=interp1(xg,double(tpl.v_grid(:))*1000,xc,'pchip',NaN);
info=struct('xQuery',xq,'xRegistered',xq,'gapLowMm',NaN(size(xq)),...
    'noGapMv',y,'gapIncrementMv',zeros(size(xq)),...
    'gapHighMm',NaN(size(xq)),'gapClippedMm',NaN(size(xq)),...
    'overshootGapMm',zeros(size(xq)),'overshootXmm',max(lo-xq,0)+max(xq-hi,0));
end

function [y,info]=eval_one(B,xGrid,gapGrid,gref,tpl,reg,gL,dg,xq,useSlope,mu,tauGap)
xq=xq(:); dg=double(dg); assert(isscalar(dg),'R5:ScalarDg');
if isfield(reg,'target_x_offset_mm') && isfinite(reg.target_x_offset_mm)
    xOffset=double(reg.target_x_offset_mm);
else
    xOffset=0;
end
xResponse=xq-xOffset;
xF=double(reg.x_scale).*(xResponse-double(reg.tau_mm));
xLo=min(double(tpl.x_grid)); xHi=max(double(tpl.x_grid)); xLowRaw=xq; xClip=min(max(xLowRaw,xLo),xHi);
% The measured low-speed template remains in its own OPR coordinate.
% target_offset belongs only to the calibrated response-surface map.
yL=interp1(double(tpl.x_grid(:)),double(tpl.v_grid(:))*1000,xClip,'pchip',NaN);
% The localization contract supplies one sensor-conditioned effective
% reference gap for each blade.  Do not re-introduce the calibration
% tilt path here: localization was performed against the scalar-gap
% surface, so applying mu(x-tau) only in the dynamic stage creates a
% static-template/model mismatch and contaminates dg.  The required
% clearance freedom remains the independent dynamic increment dg.
if useSlope && isfinite(mu) && isfinite(tauGap)
    gHL=gL+mu.*(xResponse-double(tauGap));
else
    gHL=repmat(gL,size(xq));
end
gHH=gHL+dg;
gMin=min(gapGrid); gMax=max(gapGrid); g0Clip=min(max(gHL,gMin),gMax); g1Clip=min(max(gHH,gMin),gMax);
f0=eval_surface(B,xGrid,gref,g0Clip,xF); f1=eval_surface(B,xGrid,gref,g1Clip,xF);
y=yL+double(reg.voltage_gain).*(f1-f0);
info=struct('xQuery',xq,'xRegistered',xF,'gapLowMm',gHL,'gapHighMm',gHH,...
    'noGapMv',yL,'gapIncrementMv',double(reg.voltage_gain).*(f1-f0),...
    'gapClippedMm',g1Clip,'overshootGapMm',max(gMin-gHH,0)+max(gHH-gMax,0),...
    'overshootXmm',max(xLo-xLowRaw,0)+max(xLowRaw-xHi,0));
end


function y=eval_surface(B,xGrid,gref,g,xq)
v=zeros(size(xq));
for k=1:3
    bk=interp1(xGrid(:),B(:,k),min(max(xq,min(xGrid)),max(xGrid)),'linear',NaN);
    if k==1, v=v+bk; elseif k==2, v=v+bk./g; else, v=v+bk.*log(g/gref); end
end
y=v;
end
