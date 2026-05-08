function templateLib = make_fixed_trust_template_library(gapList, xCell, yCell, ...
    gHoldout, xGridN, modelDef, trustInfo)
%make_fixed_trust_template_library  Build a template library on saved domain.

templateLib = make_response_template_library(gapList, xCell, yCell, ...
    gHoldout, xGridN, modelDef, struct('enable', false));
templateLib = restrict_gap_template_domain(templateLib, trustInfo.domain);
templateLib.trustInfo = trustInfo;
templateLib.trustInfo.selectionMode = "fixed_saved_trust_domain";
end
