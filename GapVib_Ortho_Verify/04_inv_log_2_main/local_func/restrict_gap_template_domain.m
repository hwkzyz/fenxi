function templateLib = restrict_gap_template_domain(templateLib, trustedDomain)
%restrict_gap_template_domain  Crop template grids to a validated x-domain.

trustedDomain = sort(trustedDomain(:)).';
mask = templateLib.xGrid(:) >= trustedDomain(1) & templateLib.xGrid(:) <= trustedDomain(2);
if nnz(mask) < 4
    error('Trusted domain [%.4g, %.4g] leaves too few template points.', ...
        trustedDomain(1), trustedDomain(2));
end

templateLib.fullDomain = templateLib.domain;
templateLib.fullXGrid = templateLib.xGrid;
templateLib.fullS = templateLib.S;
templateLib.fullFxMat = templateLib.FxMat;
templateLib.xGrid = templateLib.xGrid(mask);
templateLib.S = templateLib.S(:, mask);
templateLib.FxMat = templateLib.FxMat(:, mask);
templateLib.domain = [min(templateLib.xGrid), max(templateLib.xGrid)];
templateLib.trustedDomain = templateLib.domain;
end
