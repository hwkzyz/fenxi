function report = Check_R5_Migration_20250527()
%CHECK_R5_MIGRATION_20250527 Validate inputs before running the R5 front end.
cfg = Config_20250527();
report = struct('pass', true, 'files', struct(), 'messages', {{}});
required = {'responseSurface', 'lowSpeedTemplate', 'gapLibrary'};
for k = 1:numel(required)
    key = required{k};
    exists = isfile(cfg.files.(key));
    report.files.(key) = exists;
    if ~exists
        report.pass = false;
        report.messages{end+1} = sprintf('Missing %s: %s', key, cfg.files.(key));
    end
end
if report.files.responseSurface
    s = load(cfg.files.responseSurface, '-mat');
    if ~isfield(s, 'responseSurface')
        report.pass = false;
        report.messages{end+1} = 'Response MAT does not contain variable responseSurface.';
    elseif isfield(s.responseSurface, 'bladeIds') && ~ismember(cfg.case.targetBlade, s.responseSurface.bladeIds)
        report.pass = false;
        report.messages{end+1} = 'Target blade is absent from responseSurface.bladeIds.';
    end
end
if report.pass
    report.messages{end+1} = '20250527 R5 inputs passed file and target-blade checks.';
end
end
