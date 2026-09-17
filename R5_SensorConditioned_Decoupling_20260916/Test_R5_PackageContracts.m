function report = Test_R5_PackageContracts()
%TEST_R5_PACKAGECONTRACTS Contract checks for the independent R5 package.
rootDir = fileparts(mfilename('fullpath'));
restoredefaultpath;
addpath(rootDir,'-begin');
addpath(fullfile(rootDir,'core'),'-begin');
addpath(fullfile(rootDir,'adapters'),'-begin');
verify_R5_Installation('case','all');
report = struct();
report.case20250527 = check_case(Config_20250527(), false);
report.case20251222 = check_case(Config_20251222(), false);
assert(isequal(report.case20250527.searchHz,[300 1000]));
assert(isequal(report.case20251222.searchHz,[300 1000]));
assert(~report.case20250527.useStaticGapSlopeInDynamic);
assert(~report.case20251222.useStaticGapSlopeInDynamic);
assert(strcmp(report.case20250527.staticTiltCarrier,'absolute_surface_state'));
assert(strcmp(report.case20251222.staticTiltCarrier,'absolute_surface_state'));
assert(~report.case20250527.allowDynamicTilt);
assert(~report.case20251222.allowDynamicTilt);
assert(exist('R5_Evaluate_CanonicalForward','file')==2);
assert(exist('R5_Audit_StaticObservationContract','file')==2);
assert(exist('R5_Run_SensorConditionedWindowIdentification','file')==2);
fprintf('Test_R5_PackageContracts PASS.\n');
end

function r = check_case(cfg, allowDefault)
r = struct('searchHz',double(cfg.frequency.searchHz), ...
    'useStaticGapSlopeInDynamic',logical(cfg.r5.useStaticGapSlopeInDynamic), ...
    'staticTiltCarrier',char(cfg.r5.staticTiltCarrier), ...
    'allowDynamicTilt',false);
if isfield(cfg.r5,'allowDynamicTilt')
    r.allowDynamicTilt = logical(cfg.r5.allowDynamicTilt);
elseif allowDefault
    r.allowDynamicTilt = false;
end
end
