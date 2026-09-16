function report = verify_R5_Installation(varargin)
%VERIFY_R5_INSTALLATION Check the packaged R5 runtime and declared inputs.
% No legacy directory or external raw-data directory is added to the path.
p = inputParser;
addParameter(p,'case','all');
parse(p,varargin{:});
rootDir = fileparts(mfilename('fullpath'));
coreDir = fullfile(rootDir,'core');
adapterDir = fullfile(rootDir,'adapters');
addpath(coreDir,'-begin');
addpath(adapterDir,'-begin');
requested = char(p.Results.case);
if strcmpi(requested,'all')
    cases = {'20250527','20251222'};
else
    cases = {requested};
end
report = struct('packageRoot',rootDir,'cases',[],'pass',true);
for k = 1:numel(cases)
    caseId = cases{k};
    switch caseId
        case '20250527'
            cfg = Config_20250527();
            required = {fullfile(rootDir,'inputs','20250527','frozen','r5_sensor_conditioned_sidecar.mat'), ...
                cfg.files.lowSpeedTemplate, ...
                fullfile(rootDir,'inputs','20250527','foundation','Result_Step05_SingleSyncDirectTemplate_B1_S136_20250527.mat')};
        case '20251222'
            cfg = Config_20251222();
            required = {fullfile(rootDir,'inputs','20251222','frozen','R5_SensorConditionedSidecar_B1_S123.mat'), ...
                cfg.files.lowSpeedTemplate, ...
                fullfile(rootDir,'inputs','20251222','foundation','Result_Step05_NoEtaVPTopKDirectTemplate_B1_S123_20251222.mat')};
        otherwise
            error('R5:UnknownCase','Unknown packaged case: %s',caseId);
    end
    missing = required(cellfun(@(f) exist(f,'file')~=2,required));
    caseReport = struct('caseId',caseId,'requiredFiles',{required}, ...
        'missingFiles',{missing},'pass',isempty(missing), ...
        'searchHz',double(cfg.frequency.searchHz), ...
        'analysisSensors',double(cfg.case.analysisSensors));
    report.cases = [report.cases; caseReport]; %#ok<AGROW>
    if ~caseReport.pass
        report.pass = false;
        error('R5:MissingPackageInput','Missing packaged %s input(s): %s', ...
            caseId,strjoin(missing,newline));
    end
end
fprintf('R5 package verification PASS (%s).\n',strjoin(cases,', '));
end
