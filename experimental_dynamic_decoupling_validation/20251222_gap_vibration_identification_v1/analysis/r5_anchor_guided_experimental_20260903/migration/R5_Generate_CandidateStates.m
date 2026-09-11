function C = R5_Generate_CandidateStates(anchor, transitions, targetStates, opts)
%R5_GENERATE_CANDIDATESTATES Create candidate waveforms from an anchor.
% The target state is deliberately enumerated here; it is selected later by
% the high-speed profile objective, never by this migration function.

arguments
    anchor (1,1) struct
    transitions (:,1) struct
    targetStates double
    opts (1,1) struct = struct()
end

opts = apply_defaults(opts, struct('topKDonor', 3, 'distanceFloor', 1e-9, ...
    'fixedScaleMu', 0, 'fixedScaleSigma', 1, 'maxManifoldDistance', inf));

x = double(anchor.x_mm(:));
ya = double(anchor.y_mv(:));
if numel(x) ~= numel(ya)
    error('R5:AnchorSchema', 'anchor.x_mm and anchor.y_mv must have equal length.');
end

validAnchor = isfinite(x) & isfinite(ya);
if isfield(anchor, 'supportMask')
    validAnchor = validAnchor & logical(anchor.supportMask(:));
end

dist = nan(numel(transitions), 1);
for k = 1:numel(transitions)
    [yk, trSupport] = resample_transition(transitions(k), 'anchorY_mv', x);
    mask = validAnchor & trSupport & isfinite(yk);
    if nnz(mask) < 3
        continue
    end
    zq = (ya(mask) - opts.fixedScaleMu) ./ max(opts.fixedScaleSigma, eps);
    zk = (yk(mask) - opts.fixedScaleMu) ./ max(opts.fixedScaleSigma, eps);
    dist(k) = sqrt(mean((zq - zk).^2));
end
if ~any(isfinite(dist))
    error('R5:NoDonor', 'No donor transition overlaps the anchor support.');
end

targetStates = double(targetStates(:).');
C = repmat(struct('stateAxisType', '', 'stateValue', NaN, ...
    'waveform_mv', [], 'supportMask', [], 'donorIndices', [], ...
    'donorDistances', [], 'manifoldDistance', NaN, 'x_mm', x, 'status', ''), ...
    numel(targetStates), 1);

for js = 1:numel(targetStates)
    eligibleAll = find(isfinite(dist) & arrayfun(@(q) ...
        transitions(q).targetState == targetStates(js), 1:numel(transitions)).');
    [~, localOrder] = sort(dist(eligibleAll), 'ascend');
    eligible = eligibleAll(localOrder(1:min(opts.topKDonor, numel(localOrder))));
    if isempty(eligible)
        C(js).status = 'unsupported_target_state';
        continue
    end
    dd = dist(eligible);
    w = 1 ./ max(dd, opts.distanceFloor);
    w = w ./ sum(w);
    delta = zeros(size(ya));
    support = validAnchor;
    for ik = 1:numel(eligible)
        tr = transitions(eligible(ik));
        [dY, trSupport] = resample_transition(tr, 'deltaY_mv', x);
        mask = support & trSupport & isfinite(dY);
        delta(mask) = delta(mask) + w(ik) .* dY(mask);
        support = support & trSupport;
    end
    yhat = ya + delta;
    donorMean = zeros(size(ya));
    for ik = 1:numel(eligible)
        [yt, ~] = resample_transition(transitions(eligible(ik)), 'targetY_mv', x);
        donorMean = donorMean + w(ik) .* yt;
    end
    manifoldDistance = sqrt(mean((yhat(support) - donorMean(support)).^2, 'omitnan'));
    C(js).stateAxisType = transitions(eligible(1)).stateAxisType;
    C(js).stateValue = targetStates(js);
    C(js).waveform_mv = yhat;
    C(js).supportMask = support;
    C(js).donorIndices = eligible(:).';
    C(js).donorDistances = dd(:).';
    C(js).manifoldDistance = manifoldDistance;
    if nnz(support) < 3
        C(js).status = 'insufficient_support';
    elseif manifoldDistance > opts.maxManifoldDistance
        C(js).status = 'manifold_gate_failed';
    else
        C(js).status = 'candidate';
    end
end
end

function [yq, support] = resample_transition(tr, field, xq)
y = double(tr.(field)(:));
if isfield(tr, 'x_mm') && numel(tr.x_mm) == numel(y)
    xt = double(tr.x_mm(:));
    yq = interp1(xt, y, xq, 'linear', NaN);
    sourceSupport = logical(tr.supportMask(:));
    support = interp1(xt, double(sourceSupport), xq, 'nearest', 0) > 0;
else
    if numel(y) ~= numel(xq)
        error('R5:TransitionSchema', 'Transition lacks x_mm required for resampling.');
    end
    yq = y; support = logical(tr.supportMask(:));
end
end

function s = apply_defaults(s, d)
names = fieldnames(d);
for i = 1:numel(names)
    if ~isfield(s, names{i}) || isempty(s.(names{i}))
        s.(names{i}) = d.(names{i});
    end
end
end
