function A = R5_Audit_SensorRoles_20241106()
%R5_AUDIT_SENSORROLES_20241106 Confirm analysis/gap sensor separation.
cfg = Config_20241106();
A = struct();
A.analysisSensors = cfg.case.analysisSensors;
A.gapSensors = cfg.case.gapSensors;
A.analysisOnlySensors = setdiff(cfg.case.analysisSensors, cfg.case.gapSensors);
A.pass = ~isempty(A.analysisOnlySensors) && ...
    all(ismember(A.gapSensors, A.analysisSensors));
A.note = 'S2 remains in the analysis chain; S5/S7 enter gap-aware localization.';
end
