function out = R5_ProfileFit_StateVibration(window, candidates, opts)
%R5_PROFILEFIT_STATEVIBRATION Profile-fit vibration for every state candidate.
% window fields: x_mm, t_s, v_mv, and optional validMask. Candidate waveforms
% are evaluated on candidate.x_mm when supplied, otherwise window.x_mm.

arguments
    window (1,1) struct
    candidates (:,1) struct
    opts (1,1) struct = struct()
end

opts = apply_defaults(opts, struct('frequencyGridHz', linspace(300, 1000, 141), ...
    'dxBoundsMm', [-0.25 0.25], 'amplitudeBoundsMm', [0 0.50], ...
    'phaseBoundsRad', [-pi pi], 'maxIterations', 400, ...
    'minSupportFraction', 0.60, 'maxProfileSamples', 1500));

required = {'x_mm', 't_s', 'v_mv'};
for i = 1:numel(required)
    if ~isfield(window, required{i})
        error('R5:WindowSchema', 'window.%s is required.', required{i});
    end
end
x = double(window.x_mm(:));
t = double(window.t_s(:));
v = double(window.v_mv(:));
mask = isfinite(x) & isfinite(t) & isfinite(v);
if isfield(window, 'validMask')
    mask = mask & logical(window.validMask(:));
end
if nnz(mask) < 8
    error('R5:WindowSchema', 'Too few valid high-speed points for profile fitting.');
end
x = x(mask); t = t(mask); v = v(mask);
if numel(x) > opts.maxProfileSamples
    keep = round(linspace(1, numel(x), opts.maxProfileSamples));
    x = x(keep); t = t(keep); v = v(keep);
end

out = repmat(struct('stateValue', NaN, 'stateAxisType', '', ...
    'J', Inf, 'EO', NaN, 'frequencyHz', NaN, 'amplitudeMm', NaN, ...
    'phaseRad', NaN, 'dxMm', NaN, 'rmseMv', Inf, 'supportFraction', 0, ...
    'status', 'invalid'), numel(candidates), 1);

for js = 1:numel(candidates)
    c = candidates(js);
    if ~isfield(c, 'waveform_mv') || isempty(c.waveform_mv) || ...
            ~strcmpi(get_field_text(c, 'status', 'candidate'), 'candidate')
        out(js).stateValue = get_field_num(c, 'stateValue', NaN);
        out(js).stateAxisType = get_field_text(c, 'stateAxisType', '');
        out(js).status = 'candidate_unavailable';
        continue
    end
    xc = x;
    if isfield(c, 'x_mm') && numel(c.x_mm) > 1
        xc = double(c.x_mm(:));
    end
    yc = double(c.waveform_mv(:));
    if numel(xc) ~= numel(yc)
        error('R5:CandidateSchema', 'Candidate x_mm and waveform_mv mismatch.');
    end
    best = struct('J', Inf, 'frequencyHz', NaN, 'A', NaN, 'phi', NaN, 'dx', NaN);
    for frequencyHz = opts.frequencyGridHz(:).'
        w = 2*pi*frequencyHz;
        z0 = [mean(opts.amplitudeBoundsMm), 0, mean(opts.dxBoundsMm)];
        objective = @(z) profile_sse(z, frequencyHz, w, x, t, v, xc, yc, ...
            opts.amplitudeBoundsMm, opts.phaseBoundsRad, opts.dxBoundsMm);
        z = fminsearch(objective, z0, optimset('MaxIter', opts.maxIterations, ...
            'Display', 'off'));
        J = objective(z);
        if J < best.J
            best = struct('J', J, 'frequencyHz', frequencyHz, 'A', clamp(z(1), opts.amplitudeBoundsMm), ...
                'phi', wrap_phase(clamp(z(2), opts.phaseBoundsRad)), ...
                'dx', clamp(z(3), opts.dxBoundsMm));
        end
    end
    yhat = evaluate_waveform(best, x, t, xc, yc);
    support = isfinite(yhat);
    out(js).stateValue = get_field_num(c, 'stateValue', NaN);
    out(js).stateAxisType = get_field_text(c, 'stateAxisType', '');
    out(js).J = best.J;
    out(js).frequencyHz = best.frequencyHz;
    rotFreqHz = get_field_num(window, 'rotFreqHz', get_field_num(window, 'rpm_hz', NaN));
    out(js).EO = best.frequencyHz / rotFreqHz;
    out(js).amplitudeMm = best.A;
    out(js).phaseRad = best.phi;
    out(js).dxMm = best.dx;
    out(js).rmseMv = sqrt(best.J / max(nnz(support), 1));
    out(js).supportFraction = nnz(support) / numel(x);
    if out(js).supportFraction >= opts.minSupportFraction && isfinite(best.J)
        out(js).status = 'profile_fitted';
    else
        out(js).status = 'support_gate_failed';
    end
end

valid = strcmp({out.status}, 'profile_fitted');
if any(valid)
    idx = find(valid);
    [~, k] = min([out(idx).J]);
    outSelected = out(idx(k));
else
    outSelected = struct('stateValue', NaN, 'stateAxisType', '', 'J', Inf, ...
        'EO', NaN, 'frequencyHz', NaN, 'amplitudeMm', NaN, 'phaseRad', NaN, ...
        'dxMm', NaN, 'rmseMv', Inf, 'supportFraction', 0, 'status', 'no_valid_state');
end
out = struct('stateProfiles', out, 'selected', outSelected, ...
    'x_mm', x, 't_s', t, 'v_mv', v);
if isfinite(outSelected.J)
    c = candidates(find([out.stateProfiles.stateValue] == outSelected.stateValue, 1));
    xc = get_candidate_x(c, x);
    out.selected.x_mm = x;
    out.selected.t_s = t;
    out.selected.observed_mv = v;
    out.selected.predicted_mv = evaluate_waveform(struct('frequencyHz',outSelected.frequencyHz, ...
        'A',outSelected.amplitudeMm,'phi',outSelected.phaseRad,'dx',outSelected.dxMm), ...
        x, t, xc, double(c.waveform_mv(:)));
end
end

function J = profile_sse(z, ~, w, x, t, v, xc, yc, boundsA, boundsPhi, boundsDx)
z(1) = clamp(z(1), boundsA);
z(2) = clamp(z(2), boundsPhi);
z(3) = clamp(z(3), boundsDx);
u = z(1) .* sin(w .* t + z(2));
y = interp1(xc, yc, x - z(3) - u, 'pchip', NaN);
r = v - y;
J = sum(r(isfinite(r)).^2);
if ~isfinite(J) || isempty(J)
    J = realmax;
end
end

function xc = get_candidate_x(c, fallback)
if isfield(c,'x_mm') && numel(c.x_mm)>1, xc = double(c.x_mm(:)); else, xc = fallback; end
end

function y = evaluate_waveform(p, x, t, xc, yc)
y = interp1(xc, yc, x - p.dx - p.A .* sin(2*pi*p.frequencyHz.*t + p.phi), 'pchip', NaN);
end

function y = clamp(x, b)
y = min(max(x, b(1)), b(2));
end

function p = wrap_phase(p)
p = mod(p + pi, 2*pi) - pi;
end

function value = get_field_num(s, name, fallback)
if isfield(s, name), value = double(s.(name)); else, value = fallback; end
end

function value = get_field_text(s, name, fallback)
if isfield(s, name), value = char(string(s.(name))); else, value = fallback; end
end

function s = apply_defaults(s, d)
names = fieldnames(d);
for i = 1:numel(names)
    if ~isfield(s, names{i}) || isempty(s.(names{i}))
        s.(names{i}) = d.(names{i});
    end
end
end
