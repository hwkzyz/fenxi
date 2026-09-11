%% 02: One continuous static waveform surface S_b(g,xi) per blade
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath')); outDir = fullfile(thisDir,'results');
familyFile = fullfile(outDir,'StaticWaveformFamily_20251222.mat');
if ~isfile(familyFile), error('Run 01_BuildStaticWaveformFamily first.'); end
S = load(familyFile,'StaticWaveformFamily'); F = S.StaticWaveformFamily;
xi = F.xiAxisMm(:); Y = F.waveformsPeakAlignedMv; Yref = F.waveformsB2FixedPeakAlignedMv;
blades = F.bladeIds(:).'; gaps = F.nominalGapMm(:).';
g0Mm = 1.0; nX = numel(xi); nB = numel(blades); nG = numel(gaps);
coeff = nan(nX,3,nB); Yfit = nan(size(Y)); effectiveMask = false(nX,nB);
coeffReference = nan(nX,3,nB); YfitReference = nan(size(Yref)); effectiveMaskReference = false(nX,nB);
for ib = 1:nB
    envelope = max(Y(:,ib,:),[],3,'omitnan');
    effectiveMask(:,ib) = envelope >= .12 * max(envelope,[],'omitnan');
    ii = find(effectiveMask(:,ib)); if ~isempty(ii), effectiveMask(ii(1):ii(end),ib)=true; end
    A = [ones(nG,1), 1./gaps(:), log(gaps(:)./g0Mm)];
    for ix = 1:nX
        if ~effectiveMask(ix,ib), continue; end
        y = squeeze(Y(ix,ib,:)); valid = isfinite(y);
        if nnz(valid) < 4, continue; end
        c = A(valid,:) \ y(valid); coeff(ix,:,ib) = c.'; Yfit(ix,ib,:) = A*c;
        yr = squeeze(Yref(ix,ib,:)); vr = isfinite(yr);
        if nnz(vr) >= 4
            cr = A(vr,:) \ yr(vr); coeffReference(ix,:,ib) = cr.'; YfitReference(ix,ib,:) = A*cr;
            effectiveMaskReference(ix,ib) = true;
        end
    end
end
for ib = 1:nB
    ii = find(effectiveMaskReference(:,ib));
    if ~isempty(ii), effectiveMaskReference(ii(1):ii(end),ib) = true; end
end
BladeWaveformSurface = struct('method','per_blade_gap_surface_in_peak_aligned_coordinate', ...
    'basis','1, 1/g, log(g/g0)', 'g0Mm',g0Mm, 'bladeIds',blades, ...
    'nominalGapMm',gaps, 'xiAxisMm',xi, 'coeff',coeff, 'effectiveMask',effectiveMask, ...
    'observedWaveformsMv',Y, 'fittedWaveformsMv',Yfit, ...
    'coordinateReference','B2 q=0 SG/quadratic peak anchored x grid; other blades not peak aligned', ...
    'coeffReference',coeffReference, 'effectiveMaskReference',effectiveMaskReference, ...
    'observedWaveformsReferenceMv',Yref, 'fittedWaveformsReferenceMv',YfitReference, ...
    'referenceAxisMm',F.xB2OprAxisMm, ...
    'sourceFamilyFile',familyFile);
save(fullfile(outDir,'BladeWaveformSurface_20251222.mat'),'BladeWaveformSurface','-v7.3');
fprintf('Per-blade waveform surfaces saved.\n');
