function A=R5_Audit_SensorConditionedDynamic_20250527(resultFile,foundationFile,outputCsv)
cfg=Config_20250527();
if nargin<1||isempty(resultFile),resultFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_result.mat');end
if nargin<2||isempty(foundationFile),foundationFile=fullfile(cfg.paths.preparedInputs,'foundation','Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat');end
if nargin<3||isempty(outputCsv),outputCsv=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_audit.csv');end
S=load(R5_StageMatForMatlab(resultFile,'r5_dyn_audit'),'out'); R=S.out.rows;
F=load(R5_StageMatForMatlab(foundationFile,'r5_foundation_audit'),'Result'); T=F.Result.Trend;
n=numel(R); dg=reshape([R.delta_gap_mm],3,n).'; eo=[R.EO].'; Aamp=[R.amplitude_mm].'; dx=[R.dx_mm].'; rmse=[R.rmse_mv].';
base=nan(n,1); for i=1:n, k=find(T.window_id==R(i).window_id,1); if ~isempty(k), base(i)=1000*T.plain_voltage_rmse(k); end,end
sensorRows=table((1:n).',[R.window_id].',eo,Aamp,dx,dg(:,1),dg(:,2),dg(:,3),rmse,base,rmse-base,[R.hit_boundary].',... 
 'VariableNames',{'row','window_id','EO','amplitude_mm','dx_mm','dg_S1_mm','dg_S3_mm','dg_S6_mm','r5_rmse_mV','foundation_rmse_mV','delta_rmse_mV','hit_boundary'});
sensorRows.rmse_S1_mV=nan(n,1); sensorRows.rmse_S3_mV=nan(n,1); sensorRows.rmse_S6_mV=nan(n,1);
for i=1:n, q=R(i).sensor_rmse_mv; sensorRows.rmse_S1_mV(i)=q(1); sensorRows.rmse_S3_mV(i)=q(2); sensorRows.rmse_S6_mV(i)=q(3); end
writetable(sensorRows,outputCsv);
A=struct('schema','R5_SENSOR_CONDITIONED_DYNAMIC_AUDIT_V1','rows',sensorRows,...
 'pass_windows',sum(strcmp({R.status},'pass')),'total_windows',n,...
 'EO_unique',unique(eo).','dg_mean_mm',mean(dg,1),'dg_std_mm',std(dg,0,1),...
 'rmse_median_mV',median(rmse),'rmse_range_mV',[min(rmse) max(rmse)],...
 'boundary_count',sum([R.hit_boundary]),'outputCsv',outputCsv);
fprintf('R5 dynamic audit: %d/%d pass; EO=%s; RMSE median %.3f mV; boundary %d.\n',...
 A.pass_windows,n,mat2str(A.EO_unique),A.rmse_median_mV,A.boundary_count);
end
