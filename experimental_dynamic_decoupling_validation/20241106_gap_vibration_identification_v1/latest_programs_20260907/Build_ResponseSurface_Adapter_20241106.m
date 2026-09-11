function outputFile = Build_ResponseSurface_Adapter_20241106(sourceFile, outputFile)
%BUILD_RESPONSESURFACE_ADAPTER_20241106 Expose the legacy bank in R5 schema.
% The legacy bank already contains GapCalibrationBank.responseSurface; this
% adapter only unwraps it and preserves the original source file.
cfg = Config_20241106();
if nargin < 1 || isempty(sourceFile)
    sourceFile = fullfile(cfg.paths.calibrationInputs, 'gap', ...
        'GapCalibrationBank_20241106_B1toB6_S57.mat');
end
if nargin < 2 || isempty(outputFile)
    outputFile = cfg.files.responseSurface;
end
S = load(sourceFile, 'GapCalibrationBank');
assert(isfield(S, 'GapCalibrationBank') && isfield(S.GapCalibrationBank, 'responseSurface'), ...
    'R5:LegacyBankMissingResponseSurface', ...
    'Legacy gap bank does not contain GapCalibrationBank.responseSurface.');
responseSurface = S.GapCalibrationBank.responseSurface;
required = {'xGrid','trueGapMm','waveforms','bladeIds'};
for k = 1:numel(required)
    assert(isfield(responseSurface, required{k}), ...
        'R5:ResponseSurfaceFieldMissing', 'Missing responseSurface.%s.', required{k});
end
responseSurface.adapterSourceFile = sourceFile;
responseSurface.adapterSchema = 'R5_RESPONSE_SURFACE_ADAPTER_V1';
if ~exist(fileparts(outputFile), 'dir'), mkdir(fileparts(outputFile)); end
save(outputFile, 'responseSurface', '-v7.3');
fprintf('20241106 response surface adapter saved: %s\n', outputFile);
end
