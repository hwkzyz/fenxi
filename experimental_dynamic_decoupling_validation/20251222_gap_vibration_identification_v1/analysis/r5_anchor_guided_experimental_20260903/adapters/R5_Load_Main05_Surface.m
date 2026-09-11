function S = R5_Load_Main05_Surface(sourceFile)
%R5_LOAD_MAIN05_SURFACE Normalize the frozen Main05 response surface.

if nargin < 1 || isempty(sourceFile)
    error('R5:MissingSurface', 'Provide the read-only Main05 response-surface MAT file.');
end
if ~isfile(sourceFile)
    error('R5:MissingSurface', 'Response surface not found: %s', sourceFile);
end

raw = load(sourceFile);
if isfield(raw, 'responseSurface')
    R = raw.responseSurface;
elseif isfield(raw, 'CorrectedGapLibrary') && ...
        isfield(raw.CorrectedGapLibrary, 'responseSurface')
    R = raw.CorrectedGapLibrary.responseSurface;
else
    error('R5:SurfaceSchema', 'No responseSurface found in %s.', sourceFile);
end

required = {'xGrid', 'waveforms', 'bladeIds'};
for i = 1:numel(required)
    if ~isfield(R, required{i})
        error('R5:SurfaceSchema', 'responseSurface.%s is missing.', required{i});
    end
end

S = R;
S.sourceFile = sourceFile;
S.x_mm = double(R.xGrid(:));
S.blade_id = double(R.bladeIds(:).');
S.waveforms_mv = double(R.waveforms);
if isfield(R, 'trueGapMm')
    S.physical_gap_mm = double(R.trueGapMm(:).');
else
    S.physical_gap_mm = [];
end
if isfield(R, 'gTrainMm')
    S.recorded_gap_axis_mm = double(R.gTrainMm(:).');
else
    S.recorded_gap_axis_mm = [];
end
if isfield(R, 'effectiveWindow')
    S.support_mask = logical(R.effectiveWindow(:));
else
    S.support_mask = all(isfinite(S.waveforms_mv), 2);
end
if isfield(R, 'voltageUnit')
    S.voltage_unit = char(string(R.voltageUnit));
else
    S.voltage_unit = 'mV';
end

if size(S.waveforms_mv, 1) ~= numel(S.x_mm)
    error('R5:SurfaceSchema', 'xGrid and waveform first dimension disagree.');
end
if size(S.waveforms_mv, 2) ~= numel(S.blade_id)
    error('R5:SurfaceSchema', 'bladeIds and waveform second dimension disagree.');
end
end
