function templateLib = make_response_template_library(gapList, xCell, yCell, ...
    gHoldout, xGridN, modelDef, trustOptions)
%make_response_template_library  Build and configure one static template set.

if nargin < 7
    trustOptions = struct();
end

templateLib = build_gap_template_library(gapList, xCell, yCell, gHoldout, xGridN);

responseModel = string(modelDef.response_model);
switch responseModel
    case "inv_g_linear"
    case "power_law"
        templateLib.gapInterpMode = "power_law";
        templateLib.gapInterpParameter = modelDef.parameter;
        templateLib.gapInterpBasisOrder = modelDef.basis_order;
    case "exponential"
        templateLib.gapInterpMode = "exponential";
        templateLib.gapInterpParameter = modelDef.parameter;
        templateLib.gapInterpBasisOrder = modelDef.basis_order;
    case "field_basis"
        templateLib.gapInterpMode = "field_basis";
        templateLib.gapInterpBasisName = string(modelDef.basis_name);
        templateLib.gapInterpBasisOrder = modelDef.basis_order;
    otherwise
        error('Unsupported response model: %s', responseModel);
end

if isfield(trustOptions, 'enable') && trustOptions.enable
    if isfield(trustOptions, 'forcedDomain') && numel(trustOptions.forcedDomain) == 2
        trustInfo = struct();
        trustInfo.domain = sort(trustOptions.forcedDomain(:)).';
        trustInfo.selectionMode = "forced_static_library_trust_domain";
    else
        trustInfo = select_gap_trust_domain(templateLib, gapList, xCell, yCell, ...
            gHoldout, xGridN, trustOptions);
    end
    templateLib = restrict_gap_template_domain(templateLib, trustInfo.domain);
    templateLib.trustInfo = trustInfo;
end
templateLib = prepare_gap_template_library(templateLib);
end
