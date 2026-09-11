% Build a diagnostic-only Step05 result whose final bundle is the complete
% core bundle.  This isolates OPR angle mapping from phase-safe expansion.
packageRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sourceFile = fullfile(packageRoot, 'results', 'fixed_gap', ...
    'diag_opr6_unequal_fullbundle', '20250526_2500-3500_t400', ...
    'FixedGap_B1_S136_diag_opr6_unequal_fullbundle.mat');
outputDir = fullfile(packageRoot, 'results', 'fixed_gap', ...
    'diag_opr6_unequal_coreonly', '20250526_2500-3500_t400');
outputFile = fullfile(outputDir, ...
    'FixedGap_B1_S136_diag_opr6_unequal_coreonly.mat');

loaded = load(sourceFile, 'Result');
Result = loaded.Result;
for iw = 1:numel(Result.WindowResult)
    wr = Result.WindowResult(iw);
    if isempty(wr.CoreBundlePreview)
        error('Window %d has no CoreBundlePreview.', iw);
    end
    Result.WindowResult(iw).BundlePreview = wr.CoreBundlePreview;
    Result.WindowResult(iw).ExpandedBundlePreview = [];
end
Result.CreatedBy = [Result.CreatedBy, ...
    ' + diagnostic core-only bundle isolation'];
Result.RunInfo.DiagnosticBundlePolicy = 'complete_core_bundle_only';
Result.RunInfo.DiagnosticSourceFile = sourceFile;

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
save(outputFile, 'Result', '-v7.3');
fprintf('Saved diagnostic core-only foundation result:\n  %s\n', outputFile);
