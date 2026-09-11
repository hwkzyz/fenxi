function A=R5_Audit_SensorSensitivity_20250527(sidecarFile,targetBlade,outputCsv)
cfg=Config_20250527(); if nargin<1||isempty(sidecarFile),sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat');end
if nargin<2||isempty(targetBlade),targetBlade=cfg.case.targetBlade;end
if nargin<3||isempty(outputCsv),outputCsv=fullfile(cfg.paths.results,'r5_sensor_sensitivity_audit.csv');end
S=load(R5_StageMatForMatlab(sidecarFile,'sens_audit'),'SensorConditionedLibrary'); L=S.SensorConditionedLibrary;
sid=[]; z=[]; gap=[]; gain=[]; rmse=[]; med=[]; p90=[];
for i=1:numel(L.sensor)
 s=L.sensor{i}; j=find([s.state.blade_id]==targetBlade,1); if isempty(j),continue,end
 B=s.state(j).B; g=s.state(j).gap_mm; dF=-B(:,2)./g.^2+B(:,3)./g; w=isfinite(dF);
 sid(end+1,1)=s.sensor_id; z(end+1,1)=s.state(j).z; gap(end+1,1)=g; gain(end+1,1)=s.registration.voltage_gain; rmse(end+1,1)=s.state(j).rmse_mv;
 med(end+1,1)=median(abs(dF(w))); p90(end+1,1)=prctile(abs(dF(w)),90);
end
T=table(sid,z,gap,gain,rmse,med,p90,'VariableNames',{'sensor_id','latent','gap_mm','voltage_gain','anchor_rmse_mv','median_abs_dFdg','p90_abs_dFdg'}); writetable(T,outputCsv);
A=struct('schema','R5_SENSOR_SENSITIVITY_AUDIT_V1','table',T,'outputCsv',outputCsv);
disp(T);
end
