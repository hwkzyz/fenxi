function T = R5_V5_LowHigh_Closure(dynamicMapFile, windowIds, outputDir, options)
%R5_V5_LOWHIGH_CLOSURE Explicit, sidecar-only high-speed batch protocol.
arguments, dynamicMapFile (1,:) char; windowIds double; outputDir (1,:) char; options.TargetBlade = []; options.SensorId = []; end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
% Static V1/V2/V4 are a separate library-level validation and are not
% recomputed for every high-speed window.
R5_Run_All_20251222('DynamicMapFile',dynamicMapFile,'OutputDir',outputDir, ...
    'RunStaticValidation',true,'TargetBlade',options.TargetBlade,'SensorId',options.SensorId);
rows=cell(numel(windowIds),1);
for k=1:numel(windowIds)
  R=R5_Run_All_20251222('DynamicMapFile',dynamicMapFile,'WindowId',windowIds(k), ...
    'OutputDir',outputDir,'RunStaticValidation',false,'TargetBlade',options.TargetBlade,'SensorId',options.SensorId);
  p=R.profile.selected; g=R.identifiabilityGate;
  rows{k}=table(windowIds(k),p.stateValue,string(p.stateAxisType),p.EO,p.frequencyHz,p.amplitudeMm,p.rmseMv, ...
    p.supportFraction,string(g.status),string(strjoin(g.reasons,'|')), ...
    'VariableNames',{'window_id','state_value','state_axis_type','EO','frequency_hz','amplitude_mm','rmse_mv','support_fraction','decision','diagnostic_reasons'});
end
T=vertcat(rows{:}); writetable(T,fullfile(outputDir,'R5_V5_low_high_closure.csv'));
end
