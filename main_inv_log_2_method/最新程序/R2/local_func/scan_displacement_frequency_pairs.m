function scan = scan_displacement_frequency_pairs(t, u, weight, cfg)
%SCAN_DISPLACEMENT_FREQUENCY_PAIRS  Weighted full-grid harmonic scan.
%
% This is an inverse-map seed generator.  It scans the configured frequency
% grids and returns several frequency pairs, rather than committing to the
% first projected candidate.

t = t(:); u = u(:);
if nargin < 3 || isempty(weight), weight = ones(size(u)); end
weight = max(weight(:), 0);
valid = isfinite(t) & isfinite(u) & isfinite(weight) & weight > 0;
t = t(valid); u = u(valid); weight = weight(valid);
if numel(t) < 10, error('Too few valid inverse displacement samples.'); end

f1Grid = cfg.f1Grid(:).';
f2Grid = cfg.f2Grid(:).';
maxKeep = get_cfg(cfg, 'inverseMapCandidateCount', 12);
minSep = get_cfg(cfg, 'inverseMapMinFrequencySeparation', 20);

W = sqrt(weight);
uMean = sum(weight .* u) / max(sum(weight), eps);
y = W .* (u - uMean);
rows = zeros(numel(f1Grid) * numel(f2Grid), 8);
n = 0;
for i = 1:numel(f1Grid)
    for j = 1:numel(f2Grid)
        f1 = f1Grid(i); f2 = f2Grid(j);
        if abs(f2 - f1) < minSep, continue; end
        Phi = [W, W .* sin(2*pi*f1*t), W .* cos(2*pi*f1*t), ...
            W .* sin(2*pi*f2*t), W .* cos(2*pi*f2*t)];
        beta = Phi \ y;
        r = y - Phi * beta;
        n = n + 1;
        rows(n,:) = [sum(r.^2), f1, f2, beta(2), beta(3), beta(4), beta(5), beta(1)];
    end
end
rows = rows(1:n,:);
[~, order] = sort(rows(:,1), 'ascend');
rows = rows(order,:);
keep = false(size(rows,1),1);
for k = 1:size(rows,1)
    if sum(keep) >= maxKeep, break; end
    if ~any(abs(rows(keep,2)-rows(k,2)) < minSep & abs(rows(keep,3)-rows(k,3)) < minSep)
        keep(k) = true;
    end
end
rows = rows(keep,:);

scan = struct();
scan.table = rows;
scan.sse = rows(:,1);
scan.f = rows(:,2:3);
scan.beta = rows(:,4:8);
scan.u_mean = uMean;
scan.valid_count = numel(u);
scan.valid_weight_fraction = mean(valid);
end

function value = get_cfg(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
