function Phi = build_gap_response_basis(g, templateLib)
%BUILD_GAP_RESPONSE_BASIS Build the configured parametric gap basis.

g = g(:);
mode = string(templateLib.gapInterpMode);
switch mode
    case "power_law"
        p = templateLib.gapInterpParameter;
        order = get_field(templateLib, 'gapInterpBasisOrder', 1);
        if order == 1
            Phi = [ones(size(g)), g.^(-p)];
        else
            Phi = [ones(size(g)), g.^(-p), g.^(-(p + 1))];
        end
    case "exponential"
        lambda = templateLib.gapInterpParameter;
        Phi = [ones(size(g)), exp(-g ./ lambda)];
    case "field_basis"
        basisName = string(get_field(templateLib, 'gapInterpBasisName', "cap_edge_1"));
        switch basisName
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
    otherwise
        error('Gap mode %s does not use a parametric response basis.', mode);
end
end

function value = get_field(S, name, defaultValue)
if isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = defaultValue;
end
end
