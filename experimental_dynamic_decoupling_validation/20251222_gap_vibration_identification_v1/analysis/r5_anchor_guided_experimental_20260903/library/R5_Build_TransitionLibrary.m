function L = R5_Build_TransitionLibrary(S, stateAxisType, stateValues)
%R5_BUILD_TRANSITIONLIBRARY Build read-only pairwise waveform increments.
% S is the normalized output of R5_Load_Main05_Surface.  For non-B2 data,
% stateValues are labels on an effective-state axis, not physical gaps.

arguments
    S (1,1) struct
    stateAxisType (1,:) char {mustBeMember(stateAxisType, ...
        {'physical_gap_mm','effective_state'})}
    stateValues double = []
end

Y = double(S.waveforms_mv);
if ndims(Y) ~= 3
    error('R5:TransitionSchema', 'Expected waveforms as x-by-blade-by-state.');
end
nX = size(Y, 1);
nBlade = size(Y, 2);
nState = size(Y, 3);
if isempty(stateValues)
    if strcmp(stateAxisType, 'physical_gap_mm') && ~isempty(S.physical_gap_mm)
        stateValues = S.physical_gap_mm;
    else
        stateValues = 1:nState;
    end
end
stateValues = double(stateValues(:).');
if numel(stateValues) ~= nState
    error('R5:TransitionSchema', 'stateValues must match waveform state dimension.');
end

if isfield(S, 'support_mask') && numel(S.support_mask) == nX
    support = logical(S.support_mask(:));
else
    support = all(isfinite(Y), 2);
end

rows = cell(nBlade * nState * max(nState - 1, 0), 1);
ir = 0;
for ib = 1:nBlade
    for ia = 1:nState
        ya = Y(:, ib, ia);
        for it = 1:nState
            if it == ia
                continue
            end
            yt = Y(:, ib, it);
            valid = support & isfinite(ya) & isfinite(yt);
            if nnz(valid) < 3
                continue
            end
            ir = ir + 1;
            rows{ir} = struct( ...
                'sourceBladeId', S.blade_id(ib), ...
                'stateAxisType', stateAxisType, ...
                'anchorState', stateValues(ia), ...
                'targetState', stateValues(it), ...
                'x_mm', S.x_mm(:), ...
                'deltaY_mv', yt - ya, ...
                'supportMask', valid, ...
                'anchorY_mv', ya, ...
                'targetY_mv', yt);
        end
    end
end
rows = rows(1:ir);
if isempty(rows)
    error('R5:TransitionEmpty', 'No valid transition was built.');
end
L = [rows{:}].';
end
