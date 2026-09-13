function audit = R5_Audit_V1_R5_IncrementBridge_20260912(oldFoundationFile, sidecarFile, templateFile, modelBuilder, outputFile, targetBlade)
% Compare V1 and R5 gap increments at identical frozen V1 dynamic states.
% This is diagnostic only: no strain, EO lock, or parameter selection.
if nargin<5, outputFile=''; end
if nargin<6 || isempty(targetBlade)
    Z=load(R5_StageMatForMatlab(sidecarFile,'r5_sidecar'),'SensorConditionedLibrary');
    targetBlade=double(Z.SensorConditionedLibrary.sensor{1}.state(1).blade_id);
end
Q=load(oldFoundationFile,'Result'); W=Q.Result.WindowResult; nW=numel(W);
rows=repmat(struct('window_id',NaN,'sensor_id',NaN,'v1_rms_mv',NaN, ...
    'r5_rms_mv',NaN,'difference_rms_mv',NaN,'gain_ratio',NaN, ...
    'correlation',NaN,'finite_fraction',NaN),nW*3,1); krow=0;
M=modelBuilder(sidecarFile,templateFile,targetBlade);
for iw=1:nW
 B0=W(iw).bundle; fit=W(iw).modelFits.gap_only; ids=double(B0.sensorIds(:).');
 for i=1:numel(ids)
  q=B0.sensorIndex==i; x=B0.X(q)-fit.dxMm-fit.uMm(q)-fit.sensorEtaMm(i);
  c=B0.CorrectedGapLibrary.sensor(i); R=B0.responseSurface;
  dg=fit.deltaGapMm(i); xw=x-fit.deltaTauMm(i);
  [v1,~]=v1_increment(R,c,dg,xw); [~,info]=M(i).evaluate(dg,x);
  r5=double(info.gapIncrementMv(:)); r5=r5(:); v1=v1(:); ok=isfinite(v1)&isfinite(r5);
  krow=krow+1; rows(krow).window_id=W(iw).windowId; rows(krow).sensor_id=ids(i); rows(krow).finite_fraction=nnz(ok)/max(numel(ok),1);
  if nnz(ok)>=8
   a=v1(ok); b=r5(ok); rows(krow).v1_rms_mv=sqrt(mean(a.^2)); rows(krow).r5_rms_mv=sqrt(mean(b.^2)); rows(krow).difference_rms_mv=sqrt(mean((b-a).^2));
   rows(krow).gain_ratio=(a'*b)/(a'*a+eps); rows(krow).correlation=corr(a,b);
  end
 end
end
rows=rows(1:krow); audit=struct('schema','R5_V1_R5_INCREMENT_BRIDGE_AUDIT_V1','rows',rows,'source',oldFoundationFile);
if ~isempty(outputFile), save(outputFile,'audit','-v7.3'); end
end

function [d,over]=v1_increment(R,c,dg,x)
xLib=c.xScale.*(x-c.tauMm); g0=c.g0Mm+c.muGapPerXMm.*(x-c.tauMm); g1=g0+dg;
f0=surface(R,g0,xLib); f1=surface(R,g1,xLib); d=c.voltageGain.*(f1-f0);
over=max(min(R.gTrainMm)-g1,0)+max(g1-max(R.gTrainMm),0); d(over>0)=NaN;
end
function y=surface(R,g,x)
b0=interp1(R.xGrid,R.coeff(:,1),x,'linear',NaN); b1=interp1(R.xGrid,R.coeff(:,2),x,'linear',NaN); b2=interp1(R.xGrid,R.coeff(:,3),x,'linear',NaN); y=b0+b1./g+b2.*log(g/R.g0Mm);
end
