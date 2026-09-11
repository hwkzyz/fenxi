%% Build diagnostic P0-P3 transfer libraries from the identifiability results.
clear; clc;

thisDir = fileparts(mfilename('fullpath'));
caseDir = fileparts(fileparts(thisDir));
sourceFile = fullfile(caseDir, 'results', ...
    'Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136_nob_harddomain_fixed12pctroi_formal_20260830.mat');
analysisFile = fullfile(thisDir, 'outputs', 'TransferIdentifiability_20250527.mat');
outputDir = fullfile(thisDir, 'outputs', 'nested_libraries');
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

S = load(sourceFile, 'CorrectedGapLibrary');
A = load(analysisFile, 'Nested');
modelNames = unique(A.Nested.model, 'stable');
manifestRows = cell(numel(modelNames), 1);

for im = 1:numel(modelNames)
    modelName = modelNames(im);
    CorrectedGapLibrary = S.CorrectedGapLibrary;
    rows = A.Nested(A.Nested.model == modelName, :);
    for is = 1:numel(CorrectedGapLibrary.sensor)
        sid = CorrectedGapLibrary.sensor(is).sensorId;
        row = rows(rows.sensorId == sid, :);
        assert(height(row) == 1, 'Expected one %s row for CH%d.', modelName, sid);
        CorrectedGapLibrary.sensor(is).g0Mm = row.g0Mm;
        CorrectedGapLibrary.sensor(is).tauMm = row.tauMm;
        CorrectedGapLibrary.sensor(is).xScale = row.xScale;
        CorrectedGapLibrary.sensor(is).muGapPerXMm = row.muGapPerXMm;
        CorrectedGapLibrary.sensor(is).tiltAngleDeg = atan(row.muGapPerXMm) * 180 / pi;
        CorrectedGapLibrary.sensor(is).voltageGain = row.voltageGain;
        CorrectedGapLibrary.sensor(is).voltageOffsetMv = 0;
        CorrectedGapLibrary.sensor(is).lowFitRmseMv = row.rmseMv;
        CorrectedGapLibrary.sensor(is).domainFeasible = logical(row.domainFeasible);
        CorrectedGapLibrary.sensor(is).gDomainMarginMm = row.gMarginMm;
        CorrectedGapLibrary.sensor(is).xDomainMarginMm = row.xMarginMm;
        CorrectedGapLibrary.sensor(is).activeBoundary = logical(row.activeBoundary);
    end
    CorrectedGapLibrary.method = sprintf('diagnostic_nested_transfer_%s', modelName);
    CorrectedGapLibrary.description = sprintf(['Diagnostic nested platform-transfer library %s. ' ...
        'Built from the same baseline-removed 12%%%% ROI and calibration domain as the formal no-b library.'], modelName);
    CorrectedGapLibrary.transferModel = char(modelName);
    CorrectedGapLibrary.sourceFormalLibrary = sourceFile;
    CorrectedGapLibrary.sourceIdentifiabilityAnalysis = analysisFile;
    outFile = fullfile(outputDir, sprintf('GapLibrary_20250527_%s.mat', modelName));
    save(outFile, 'CorrectedGapLibrary', '-v7');
    manifestRows{im} = table(modelName, string(outFile), ...
        mean(rows.rmseMv), all(rows.domainFeasible), ...
        'VariableNames', {'model','libraryFile','meanLowFitRmseMv','allDomainFeasible'});
    fprintf('%s -> %s\n', modelName, outFile);
end

Manifest = vertcat(manifestRows{:});
writetable(Manifest, fullfile(outputDir, 'NestedTransferLibraryManifest_20250527.csv'));
disp(Manifest);

