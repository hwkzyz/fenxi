function C = R5_Validate_InputContract(cfg, surface, template, dynamicFile)
%R5_VALIDATE_INPUTCONTRACT Check units, dimensions and source traceability.
C = struct('status','pass','messages',{{}},'sourceFiles',{{}},'surface',struct(), 'template',struct());
C.sourceFiles = {cfg.inputs.responseSurface, cfg.inputs.lowSpeedTemplate};
if ~isempty(dynamicFile), C.sourceFiles{end+1} = dynamicFile; end
assert(isfield(surface,'x_mm') && isfield(surface,'waveforms_mv'), 'R5:Contract', 'Main05 surface lacks x/waveform.');
Y = double(surface.waveforms_mv);
assert(ndims(Y)==3, 'R5:Contract', 'Main05 waveforms must be x-by-blade-by-state.');
assert(numel(surface.x_mm)==size(Y,1), 'R5:Contract', 'x grid and waveform dimensions disagree.');
assert(all(isfinite(surface.x_mm)), 'R5:Contract', 'Main05 x grid contains non-finite values.');
assert(startsWith(lower(strtrim(surface.voltage_unit)),'mv'), 'R5:Contract', 'Waveform unit must be mV.');
C.surface.nX = size(Y,1); C.surface.nBlade = size(Y,2); C.surface.nState = size(Y,3);
C.surface.supportFraction = nnz(surface.support_mask)/numel(surface.support_mask);
assert(C.surface.supportFraction >= 0.5, 'R5:Contract', 'Main05 support fraction is too small.');
if ~isempty(template)
    C.template.hasTemplate = isfield(template,'Window') || isfield(template,'Sensor') || isfield(template,'sensor');
else
    C.template.hasTemplate = false;
end
if ~isempty(dynamicFile), assert(isfile(dynamicFile), 'R5:Contract', 'DynamicMap file is not readable.'); end
end
