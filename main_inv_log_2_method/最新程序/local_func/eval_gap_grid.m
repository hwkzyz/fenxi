function [curveGrid, dGrid] = eval_gap_grid(templateLib, g)
%eval_gap_grid  Interpolate template and derivative grids along gap.

if isfield(templateLib, 'gapInterpMode')
    [gSort, order] = sort(templateLib.gapTrain(:), 'ascend');
    switch string(templateLib.gapInterpMode)
        case "g_linear"
            curveGrid = interp1(gSort, templateLib.S(order, :), g, 'linear', 'extrap');
            dGrid = interp1(gSort, templateLib.FxMat(order, :), g, 'linear', 'extrap');
            return;
        case "g_pchip"
            curveGrid = interp1(gSort, templateLib.S(order, :), g, 'pchip', 'extrap');
            dGrid = interp1(gSort, templateLib.FxMat(order, :), g, 'pchip', 'extrap');
            return;
        case "power_law"
            p = templateLib.gapInterpParameter;
            basisOrder = get_template_field(templateLib, 'gapInterpBasisOrder', 1);
            Phi = build_response_basis(gSort, "power_law", p, basisOrder);
            phiQuery = build_response_basis(g, "power_law", p, basisOrder);
            curveGrid = phiQuery * (Phi \ templateLib.S(order, :));
            dGrid = phiQuery * (Phi \ templateLib.FxMat(order, :));
            return;
        case "exponential"
            lambda = templateLib.gapInterpParameter;
            Phi = build_response_basis(gSort, "exponential", lambda, 1);
            phiQuery = build_response_basis(g, "exponential", lambda, 1);
            curveGrid = phiQuery * (Phi \ templateLib.S(order, :));
            dGrid = phiQuery * (Phi \ templateLib.FxMat(order, :));
            return;
        case "field_basis"
            basisName = get_template_field(templateLib, 'gapInterpBasisName', "cap_edge_1");
            Phi = build_field_basis(gSort, basisName);
            phiQuery = build_field_basis(g, basisName);
            curveGrid = phiQuery * (Phi \ templateLib.S(order, :));
            dGrid = phiQuery * (Phi \ templateLib.FxMat(order, :));
            return;
    end
end

zTrain = 1 ./ templateLib.gapTrain(:);
[zSort, order] = sort(zTrain, 'ascend');
zq = 1 / max(g, 1e-6);
curveGrid = interp1(zSort, templateLib.S(order, :), zq, 'linear', 'extrap');
dGrid = interp1(zSort, templateLib.FxMat(order, :), zq, 'linear', 'extrap');
end

function Phi = build_response_basis(g, modelName, parameter, basisOrder)
g = g(:);
switch string(modelName)
    case "power_law"
        if basisOrder == 1
            Phi = [ones(size(g)), g.^(-parameter)];
        else
            Phi = [ones(size(g)), g.^(-parameter), g.^(-(parameter + 1))];
        end
    case "exponential"
        Phi = [ones(size(g)), exp(-g ./ parameter)];
end
end

function Phi = build_field_basis(g, basisName)
g = g(:);
switch string(basisName)
    case "cap_edge_1"
        Phi = [ones(size(g)), 1 ./ g, log(g), 1 ./ (g .^ 2)];
    case "cap_edge_2"
        Phi = [ones(size(g)), 1 ./ g, 1 ./ (g .^ 2), log(g)];
    case "cap_edge_3"
        Phi = [ones(size(g)), 1 ./ g, log(g), g .* log(g)];
    case "inv_poly_3"
        Phi = [ones(size(g)), 1 ./ g, 1 ./ (g .^ 2), 1 ./ (g .^ 3)];
    case "inv_log_2"
        Phi = [ones(size(g)), 1 ./ g, log(g)];
    otherwise
        error('Unknown field-inspired basis: %s', basisName);
end
end

function val = get_template_field(S, name, defaultVal)
if isfield(S, name) && ~isempty(S.(name))
    val = S.(name);
else
    val = defaultVal;
end
end
