function report = Check_R5_Migration_20241106()
%CHECK_R5_MIGRATION_20241106 Report the data contract and migration blocker.
cfg = Config_20241106();
report = struct('pass', true, 'files', struct(), 'messages', {{}});
keys = {'responseSurface','lowSpeedTemplate','gapLibrary'};
for k = 1:numel(keys)
    key = keys{k}; exists = isfile(cfg.files.(key)); report.files.(key) = exists;
    if ~exists
        report.pass = false;
        report.messages{end+1} = sprintf('Missing %s: %s', key, cfg.files.(key));
    end
end
if ~report.files.responseSurface
    report.messages{end+1} = ['Run Build_ResponseSurface_Adapter_20241106 using ', ...
        'GapCalibrationBank_20241106_B1toB6_S57.mat; do not reconstruct waveforms from the CSV.'];
end
end
