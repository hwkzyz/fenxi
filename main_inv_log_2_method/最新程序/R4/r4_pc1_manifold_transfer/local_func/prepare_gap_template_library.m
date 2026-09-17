function templateLib = prepare_gap_template_library(templateLib)
%PREPARE_GAP_TEMPLATE_LIBRARY Cache fixed parametric-gap coefficients.

if ~isfield(templateLib, 'gapInterpMode')
    return;
end
mode = string(templateLib.gapInterpMode);
if ~ismember(mode, ["power_law", "exponential", "field_basis"])
    return;
end

[gSort, order] = sort(templateLib.gapTrain(:), 'ascend');
Phi = build_gap_response_basis(gSort, templateLib);
templateLib.gapBasisCurveCoefficients = Phi \ templateLib.S(order, :);
templateLib.gapBasisDerivativeCoefficients = Phi \ templateLib.FxMat(order, :);
end
