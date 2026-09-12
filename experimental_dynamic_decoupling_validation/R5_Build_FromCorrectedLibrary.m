function models = R5_Build_FromCorrectedLibrary(sidecarFile, templateFile, targetBlade)
% Common R5 model builder for converted Main10 response libraries.
S = load(R5_StageMatForMatlab(sidecarFile,'r5_sidecar'),'SensorConditionedLibrary');
load(R5_StageMatForMatlab(templateFile,'r5_template'),'Template');
L=S.SensorConditionedLibrary; assert(strcmpi(L.latent_mode,'sensor_conditioned'),'R5:LatentMode');
if isfield(Template,'Sensor'), arr=Template.Sensor; elseif isfield(Template,'SensorBlade'), arr=Template.SensorBlade; else, error('R5:TemplateSensorsMissing'); end
models=struct('sensorId',{},'bladeId',{},'gapReferenceMm',{},'latent',{},'isGapSensor',{},'coord',{},'templateDomainMm',{},'responseDomainMm',{},'gapDomainMm',{},'templateBaseline',{},'templateVoltageUnit',{},'templateEvaluate',{},'evaluate',{});
for i=1:numel(L.sensor)
 sc=L.sensor{i}; st=sc.state(1); reg=sc.registration; sid=double(sc.sensor_id); ib=find([arr.sensor_id]==sid,1); assert(~isempty(ib),'R5:TemplateMissing'); tpl=arr(ib); C=double(sc.coeff_by_state{1});
 coord=struct('templateOffsetMm',0,'responseOffsetMm',getf(reg,'target_x_offset_mm',0),'responseTauMm',getf(reg,'tau_mm',0),'responseScale',getf(reg,'x_scale',1),'responseMuGapPerXMm',getf(reg,'mu_gap_per_x_mm',0));
 models(end+1)=struct('sensorId',sid,'bladeId',double(targetBlade),'gapReferenceMm',double(st.gap_mm),'latent',double(st.z),'isGapSensor',true,'coord',coord,'templateDomainMm',[min(tpl.x_grid) max(tpl.x_grid)],'responseDomainMm',[min(L.x_mm) max(L.x_mm)],'gapDomainMm',[min(L.gap_mm) max(L.gap_mm)],'templateBaseline',double(tpl.baseline),'templateVoltageUnit','V','templateEvaluate',@(xq)et(tpl,xq),'evaluate',@(dg,xq)ev(C,L.x_mm,L.gap_mm,L.gref_mm,tpl,reg,st.gap_mm,dg,xq));
end
models=models(:).';
end
function v=getf(s,n,d), v=d; if isfield(s,n)&&isfinite(s.(n)), v=double(s.(n)); end, end
function y=et(t,x), y=interp1(double(t.x_grid(:)),double(t.v_grid(:))*1000,double(x(:)),'pchip',NaN); end
function [y,info]=ev(C,X,G,gr,t,reg,gL,dg,x)
x=double(x(:)); tau=getf(reg,'tau_mm',0); xs=getf(reg,'x_scale',1); mu=getf(reg,'mu_gap_per_x_mm',0); gain=getf(reg,'voltage_gain',1); xf=xs.*(x-tau); g0=double(gL)+mu.*(x-tau); g1=g0+double(dg); f0=es(C,X,gr,g0,xf); f1=es(C,X,gr,g1,xf); y=et(t,x)+gain.*(f1-f0); valid=x>=min(t.x_grid)&x<=max(t.x_grid)&xf>=min(X)&xf<=max(X)&g0>=min(G)&g0<=max(G)&g1>=min(G)&g1<=max(G); y(~valid)=NaN; info=struct('xQuery',x,'xRegistered',xf,'gapLowMm',g0,'gapHighMm',g1,'noGapMv',et(t,x),'gapIncrementMv',gain.*(f1-f0),'gapClippedMm',g1,'overshootGapMm',max(min(G)-g1,0)+max(g1-max(G),0),'overshootXmm',max(min(t.x_grid)-x,0)+max(x-max(t.x_grid),0),'overshootResponseXmm',max(min(X)-xf,0)+max(xf-max(X),0));
end
function y=es(C,X,gr,g,x), b0=interp1(X(:),C(:,1),x,'linear',NaN); b1=interp1(X(:),C(:,2),x,'linear',NaN); b2=interp1(X(:),C(:,3),x,'linear',NaN); y=b0+b1./g+b2.*log(g/gr); end
