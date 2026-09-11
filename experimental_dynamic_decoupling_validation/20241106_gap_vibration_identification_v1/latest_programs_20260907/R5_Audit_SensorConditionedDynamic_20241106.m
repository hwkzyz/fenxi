function A = R5_Audit_SensorConditionedDynamic_20241106(resultFile,foundationFile,outputCsv)
cfg=Config_20241106();
if nargin<1||isempty(resultFile),resultFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_result.mat');end
if nargin<2||isempty(foundationFile),foundationFile='';end
if nargin<3||isempty(outputCsv),outputCsv=fullfile(cfg.paths.results,'r5_sensor_conditioned_dynamic_audit.csv');end
S=load(R5_StageMatForMatlab(resultFile,'r5_dyn_audit'),'out'); O=S.out; R=O.rows;
n=numel(R); allSensorIds=double(O.sensorIds(:).'); nS=numel(allSensorIds);
dg=reshape([R.delta_gap_mm],nS,n).';
eo=[R.EO].'; amp=[R.amplitude_mm].'; dx=[R.dx_mm].'; rmse=[R.rmse_mv].';
V=table((1:n).',[R.window_id].',eo,amp,dx,rmse,[R.hit_boundary].', ...
    'VariableNames',{'row','window_id','EO','amplitude_mm','dx_mm','r5_rmse_mV','hit_boundary'});
for sid=cfg.case.gapSensors(:).'
    k=find(allSensorIds==sid,1);
    assert(~isempty(k),'R5:AuditSensorRole','Gap sensor S%d is missing from the result.',sid);
    V.(sprintf('dg_S%d_mm',sid))=dg(:,k);
    q=nan(n,1); for i=1:n, s=R(i).sensor_rmse_mv; if numel(s)>=k,q(i)=s(k);end,end
    V.(sprintf('rmse_S%d_mV',sid))=q;
end
if ~isempty(foundationFile)
    F=load(R5_StageMatForMatlab(foundationFile,'r5_foundation_audit'),'Result');
    if isfield(F.Result,'Trend')
        T=F.Result.Trend; base=nan(n,1);
        for i=1:n, k=find(T.window_id==R(i).window_id,1); if ~isempty(k),base(i)=1000*T.plain_voltage_rmse(k);end,end
        V.foundation_rmse_mV=base; V.delta_rmse_mV=rmse-base;
    end
end
writetable(V,outputCsv);
A=struct('schema','R5_SENSOR_CONDITIONED_DYNAMIC_AUDIT_V1','rows',V, ...
    'sensorIds',allSensorIds,'gapSensorIds',double(cfg.case.gapSensors(:).'), ...
    'pass_windows',sum(strcmp({R.status},'pass')),'total_windows',n, ...
    'EO_unique',unique(eo).','dg_mean_mm',mean(dg,1),'dg_std_mm',std(dg,0,1), ...
    'rmse_median_mV',median(rmse),'rmse_range_mV',[min(rmse) max(rmse)], ...
    'boundary_count',sum([R.hit_boundary]),'outputCsv',outputCsv);
fprintf('20241106 R5 dynamic audit: %d/%d pass; EO=%s; RMSE median %.3f mV; boundary %d.\n', ...
    A.pass_windows,n,mat2str(A.EO_unique),A.rmse_median_mV,A.boundary_count);
end
