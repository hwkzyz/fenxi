function [safeMap, info] = apply_query_safe_window_mask(highMap, templateLib, cfg)
%APPLY_QUERY_SAFE_WINDOW_MASK Freeze one safe sample set for all fit stages.
% The final query is x - dx - u(t) - eta. Remove nominal x samples that
% could leave the trusted template domain under the configured bounds.

if nargin < 3 || isempty(cfg), cfg = struct(); end
if isfield(highMap, 'query_safe_applied') && highMap.query_safe_applied
    safeMap = highMap;
    info = get_query_safe_info(highMap);
    return;
end

x = highMap.x_v(:);
n = numel(x);
xDomain = resolve_domain(templateLib);
ampBound = get_cfg(cfg, 'querySafeAmplitudeBoundMm', []);
if isempty(ampBound), ampBound = get_cfg(cfg, 'finalAmplitudeUpperBoundMm', []); end
if isempty(ampBound), ampBound = get_cfg(cfg, 'nonlinearAmpUpperBound', 0.5); end
dxBound = get_cfg(cfg, 'querySafeDxBoundMm', []);
if isempty(dxBound), dxBound = get_cfg(cfg, 'finalProjectedJointDxHalfWidth', []); end
if isempty(dxBound), dxBound = get_cfg(cfg, 'projectedJointDxHalfWidth', 0.20); end
etaBound = get_cfg(cfg, 'querySafeEtaBoundMm', 0);
margin = get_cfg(cfg, 'querySafeDomainMarginMm', 0.02);
guard = max(0, abs(ampBound) + abs(dxBound) + abs(etaBound) + abs(margin));
safeDomain = [xDomain(1) + guard, xDomain(2) - guard];
if safeDomain(2) <= safeDomain(1)
    error('invlog2:QuerySafeDomainCollapsed', ...
        'Query-safe guard %.6g mm collapses template domain [%.6g, %.6g].', ...
        guard, xDomain(1), xDomain(2));
end

mask = isfinite(x) & x >= safeDomain(1) & x <= safeDomain(2);
safeMap = highMap;
names = fieldnames(highMap);
for i = 1:numel(names)
    name = names{i};
    value = highMap.(name);
    if isnumeric(value) || islogical(value) || isstring(value) || iscell(value)
        if numel(value) == n && ~isscalar(value)
            safeMap.(name) = value(mask);
        end
    end
end
safeMap.query_safe_applied = true;
safeMap.query_safe_mask_original = mask;
safeMap.query_safe_guard_mm = guard;
safeMap.query_safe_domain = safeDomain;
safeMap.query_safe_source = 'fixed_parameter_bounds';
safeMap.query_safe_original_count = n;
safeMap.query_safe_count = nnz(mask);
safeMap.query_safe_removed_count = nnz(~mask);
info = get_query_safe_info(safeMap);
end

function domain = resolve_domain(lib)
if isfield(lib, 'domain') && numel(lib.domain) >= 2 && all(isfinite(lib.domain(1:2)))
    domain = double(lib.domain(1:2));
elseif isfield(lib, 'xGrid')
    domain = [min(lib.xGrid(:)), max(lib.xGrid(:))];
else
    error('invlog2:MissingTemplateDomain', 'Template library has no usable x-domain.');
end
end

function value = get_cfg(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name)), value = cfg.(name); else, value = defaultValue; end
end

function info = get_query_safe_info(S)
info = struct('enabled', true, 'guard_mm', NaN, 'domain', [NaN NaN], ...
    'original_count', NaN, 'safe_count', NaN, 'removed_count', NaN, 'source', '');
if isfield(S, 'query_safe_guard_mm'), info.guard_mm = S.query_safe_guard_mm; end
if isfield(S, 'query_safe_domain'), info.domain = S.query_safe_domain; end
if isfield(S, 'query_safe_original_count'), info.original_count = S.query_safe_original_count; end
if isfield(S, 'query_safe_count'), info.safe_count = S.query_safe_count; end
if isfield(S, 'query_safe_removed_count'), info.removed_count = S.query_safe_removed_count; end
if isfield(S, 'query_safe_source'), info.source = S.query_safe_source; end
end
