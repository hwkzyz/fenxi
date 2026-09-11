function models = R5_Build_SensorConditionedDynamicModels_20241106(sidecarFile, templateFile, targetBlade)
% Build the formal R5 dynamic interface: only dg_s is window-varying.
cfg = Config_20241106();
if nargin < 1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin < 2 || isempty(templateFile), templateFile=cfg.files.lowSpeedTemplate; end
if nargin < 3 || isempty(targetBlade), targetBlade=cfg.case.targetBlade; end
S=load(R5_StageMatForMatlab(sidecarFile,'r5_sidecar'),'SensorConditionedLibrary'); load(R5_StageMatForMatlab(templateFile,'r5_template'),'Template');
L=S.SensorConditionedLibrary; assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
models=struct('sensorId',{},'bladeId',{},'gapReferenceMm',{},'latent',{},'evaluate',{});
for i=1:numel(L.sensor)
    sc=L.sensor{i}; j=find([sc.state.blade_id]==targetBlade,1);
    if isempty(j), continue; end
    sid=double(sc.sensor_id); ib=find([Template.SensorBlade.sensor_id]==sid & ...
        [Template.SensorBlade.blade_id]==targetBlade,1);
    assert(~isempty(ib),'R5:TemplateMissing','Missing low-speed template for B%d/S%d.',targetBlade,sid);
    tpl=Template.SensorBlade(ib); reg=sc.registration; st=sc.state(j);
    if isfield(st,'target_x_offset_mm'), reg.target_x_offset_mm=double(st.target_x_offset_mm); end
    m=struct('sensorId',sid,'bladeId',targetBlade,'gapReferenceMm',st.gap_mm,...
        'latent',st.z,'evaluate',@(dg,xq) eval_one(st.B,L.x_mm,L.gap_mm,L.gref_mm,...
        tpl,reg,st.gap_mm,dg,xq,cfg.r5.useStaticGapSlopeInDynamic,st.mu_gap_per_x_mm,st.gap_tau_mm));
    models(end+1)=m; %#ok<AGROW>
end
assert(~isempty(models),'R5:NoDynamicSensors','No target-blade sensor states found.');
end

function [y,info]=eval_one(B,xGrid,gapGrid,gref,tpl,reg,gL,dg,xq,useSlope,mu,tauGap)
xq=xq(:); dg=double(dg); assert(isscalar(dg),'R5:ScalarDg');
xOffset=0;
xF=double(reg.x_scale).*(xq-double(reg.tau_mm));
xLo=min(double(tpl.x_grid)); xHi=max(double(tpl.x_grid)); xLowRaw=xq-xOffset; xClip=min(max(xLowRaw,xLo),xHi);
yL=interp1(double(tpl.x_grid(:)),(double(tpl.v_grid(:))-double(tpl.baseline))*1000,xClip,'pchip',NaN);
% The localization contract supplies one sensor-conditioned effective
% reference gap for each blade.  Do not re-introduce the calibration
% tilt path here: localization was performed against the scalar-gap
% surface, so applying mu(x-tau) only in the dynamic stage creates a
% static-template/model mismatch and contaminates dg.  The required
% clearance freedom remains the independent dynamic increment dg.
if useSlope && isfinite(mu) && isfinite(tauGap)
    gHL=gL+mu.*(xq-double(tauGap));
else
    gHL=repmat(gL,size(xq));
end
gHH=gHL+dg;
gMin=min(gapGrid); gMax=max(gapGrid); g0Clip=min(max(gHL,gMin),gMax); g1Clip=min(max(gHH,gMin),gMax);
f0=eval_surface(B,xGrid,gref,g0Clip,xF); f1=eval_surface(B,xGrid,gref,g1Clip,xF);
y=yL+double(reg.voltage_gain).*(f1-f0);
info=struct('xQuery',xq,'xRegistered',xF,'gapLowMm',gHL,'gapHighMm',gHH,...
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
