function report = R5_CheckPackageDependencies()
%R5_CHECKPACKAGEDEPENDENCIES Verify that formal package code is self-contained.
rootDir = fileparts(mfilename('fullpath'));
codeDirs = {fullfile(rootDir,'core'),fullfile(rootDir,'adapters'),rootDir};
files = {};
for k = 1:numel(codeDirs)
    d = dir(fullfile(codeDirs{k},'*.m'));
    files = [files; fullfile({d.folder},{d.name})']; %#ok<AGROW>
end
files = unique(files);
legacyTokens = {'latest_programs','Main06','Main10','explicit_profile','g_eff = g0 + mu'}; % checker vocabulary
externalTokens = {'E:\\试验数据','pwd','cd('};
violations = struct('file',{},'legacyTokens',{},'externalTokens',{});
for k = 1:numel(files)
    [~,fileName,fileExt] = fileparts(files{k});
    if strcmpi([fileName fileExt],'R5_CheckPackageDependencies.m')
        continue
    end
    txt = fileread(files{k});
    hitLegacy = legacyTokens(cellfun(@(s) contains(txt,s),legacyTokens));
    hitExternal = externalTokens(cellfun(@(s) contains(txt,s),externalTokens));
    if ~isempty(hitLegacy) || ~isempty(hitExternal)
        violations(end+1) = struct('file',files{k}, ...
            'legacyTokens',{hitLegacy},'externalTokens',{hitExternal}); %#ok<AGROW>
    end
end
report = struct('packageRoot',rootDir,'matlabFiles',{files}, ...
    'violations',{violations},'pass',isempty(violations));
if ~isempty(violations)
    for k = 1:numel(violations)
        fprintf('VIOLATION: %s\n',violations(k).file);
        if ~isempty(violations(k).legacyTokens)
            fprintf('  legacy: %s\n',strjoin(violations(k).legacyTokens,', '));
        end
        if ~isempty(violations(k).externalTokens)
            fprintf('  external: %s\n',strjoin(violations(k).externalTokens,', '));
        end
    end
end
assert(report.pass,'R5:PackageDependencyViolation', ...
    'Package dependency audit found forbidden references.');
fprintf('R5 package dependency audit PASS (%d MATLAB files).\n',numel(files));
end
