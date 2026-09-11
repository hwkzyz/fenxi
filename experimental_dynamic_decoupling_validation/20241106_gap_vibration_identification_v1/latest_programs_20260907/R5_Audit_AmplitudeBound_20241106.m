function summary = R5_Audit_AmplitudeBound_20241106(foundationFile, outputCsv)
%R5_AUDIT_AMPLITUDEBOUND_20241106 Separate a physical amplitude bound from model fit.
% This is a diagnostic rerun: it does not replace the formal 0.50-mm result.
rootDir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(rootDir);
cfg=Config_20241106();
if nargin<1 || isempty(foundationFile)
    error('R5:MissingFoundation','Pass the archived or current Foundation Result file.');
end
if nargin<2 || isempty(outputCsv)
    outputCsv=fullfile(cfg.paths.results,'r5_amplitude_bound_audit.csv');
end
sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat');
templateFile=fullfile(cfg.paths.results,'r5_low_speed_template_bank_adapter.mat');
bounds=[0.50 0.80 1.00];
rows=repmat(struct('amplitude_limit_mm',NaN,'window_count',0,'boundary_count',0,...
    'median_amplitude_mm',NaN,'max_amplitude_mm',NaN,'median_rmse_mv',NaN,...
    'median_foundation_rmse_mv',NaN,'median_delta_rmse_mv',NaN),numel(bounds),1);
for k=1:numel(bounds)
    cfg.r5.amplitudeLimitMm=bounds(k);
    outFile=fullfile(cfg.paths.results,sprintf('r5_dynamic_ampbound_%0.2f.mat',bounds(k)));
    out=R5_Run_SensorConditionedWindowIdentification(cfg,@R5_Build_SensorConditionedDynamicModels_20241106,...
        sidecarFile,templateFile,foundationFile,outFile);
    R=out.rows; valid=strcmp({R.status},'pass');
    rows(k).amplitude_limit_mm=bounds(k); rows(k).window_count=sum(valid);
    rows(k).boundary_count=sum([R.hit_boundary]);
    rows(k).median_amplitude_mm=median([R(valid).amplitude_mm]);
    rows(k).max_amplitude_mm=max([R(valid).amplitude_mm]);
    rows(k).median_rmse_mv=median([R(valid).rmse_mv]);
    F=load(R5_StageMatForMatlab(foundationFile,'r5_amp_audit'),'Result'); T=F.Result.Trend;
    fr=double([T.plain_voltage_rmse]); fr=fr(isfinite(fr));
    if ~isempty(fr) && median(abs(fr))<1, fr=1000*fr; end
    rows(k).median_foundation_rmse_mv=median(fr);
    rows(k).median_delta_rmse_mv=rows(k).median_rmse_mv-rows(k).median_foundation_rmse_mv;
end
summary=struct2table(rows);
if ~exist(fileparts(outputCsv),'dir'), mkdir(fileparts(outputCsv)); end
writetable(summary,outputCsv);
end
