function low = aggregate_low_speed_template(dataLow, alpha_k, R_tip, domain, xGrid)
%AGGREGATE_LOW_SPEED_TEMPLATE  Map and aggregate low-speed turns in x-space.

mapped = map_highspeed_to_space(dataLow, alpha_k, R_tip, domain, 0);
xGrid = xGrid(:);
keys = [mapped.rev_v(:), mapped.S_v(:)];
if isempty(keys)
    error('No low-speed samples fall inside the requested spatial domain.');
end
groupKeys = unique(keys, 'rows', 'stable');
Y = nan(size(groupKeys, 1), numel(xGrid));
for ig = 1:size(groupKeys, 1)
    idx = mapped.rev_v == groupKeys(ig, 1) & mapped.S_v == groupKeys(ig, 2);
    x = mapped.x_v(idx);
    v = mapped.V_a(idx);
    [x, order] = sort(x(:));
    v = v(order);
    [x, keep] = unique(x, 'stable');
    v = v(keep);
    if numel(x) >= 2
        inside = xGrid >= x(1) & xGrid <= x(end);
        Y(ig, inside) = interp1(x, v, xGrid(inside), 'pchip');
    end
end

templateLow = mean(Y, 1, 'omitnan').';
stdLow = std(Y, 0, 1, 'omitnan').';
countLow = sum(isfinite(Y), 1).';
valid = isfinite(templateLow) & countLow > 0;
if nnz(valid) < 0.9 * numel(xGrid)
    error('Insufficient low-speed spatial coverage for template aggregation.');
end

low = struct();
low.xGrid = xGrid;
low.templateLow = fillmissing(templateLow, 'linear', 'EndValues', 'nearest');
low.stdLow = fillmissing(stdLow, 'nearest');
low.countLow = countLow;
low.perTurn = Y;
low.mapped = mapped;
low.num_turns = size(Y, 1);
end
