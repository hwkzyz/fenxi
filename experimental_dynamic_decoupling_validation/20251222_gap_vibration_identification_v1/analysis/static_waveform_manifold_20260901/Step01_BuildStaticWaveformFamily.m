%% 01: Static waveform family in B2-anchored and diagnostic coordinates
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(thisDir));
addpath(rootDir);
outDir = fullfile(thisDir, 'results');
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
sourceFile = fullfile(rootDir, 'results', 'prepared', 'calibration', ...
    'Step05_Response_Surface_20251222_RefBladeAnchor.mat');
if ~isfile(sourceFile), error('Run Main05_RefBladeAnchor_20251222 first.'); end

S = load(sourceFile, 'responseSurface');
R = S.responseSurface;
x = R.xGrid(:); Yx = R.waveforms; blades = R.bladeIds(:).'; gaps = R.trueGapMm(:).';
nX = numel(x); nB = numel(blades); nG = numel(gaps);
xi = x; Yxi = nan(size(Yx)); peakXiMm = nan(nB, nG); peakVoltageMv = nan(nB, nG);
peakGridMm = nan(nB, nG);

for ib = 1:nB
    for ig = 1:nG
        y = Yx(:,ib,ig);
        valid = isfinite(y);
        if nnz(valid) < 11, continue; end
        xv = x(valid); yv = y(valid);
        ys = sgolayfilt(yv, 3, 11);
        [~, ii] = max(ys);
        lo = max(1, ii-2); hi = min(numel(xv), ii+2);
        xf = xv(lo:hi); yf = ys(lo:hi);
        xpk = xv(ii); vpk = ys(ii);
        if numel(xf) >= 3
            p = polyfit(xf, yf, 2);
            if p(1) < 0
                xq = -p(2)/(2*p(1));
                if xq >= xf(1) && xq <= xf(end)
                    xpk = xq; vpk = polyval(p, xq);
                end
            end
        end
        peakGridMm(ib,ig) = xv(ii);
        peakXiMm(ib,ig) = xpk;
        peakVoltageMv(ib,ig) = vpk;
        Yxi(:,ib,ig) = interp1(xv - xpk, yv, xi, 'pchip', NaN);
    end
end

% B2 is the q=0 reference.  Use its SG/quadratic peak as the per-state
% origin, then apply the same small coordinate correction to every blade in
% that state.  Other blades are never independently peak-aligned.
b2 = find(blades == 2, 1);
Yb2Fixed = nan(size(Yx));
for ig = 1:nG
    d = peakXiMm(b2,ig);
    Yb2Fixed(:,:,ig) = interp1(x, Yx(:,:,ig), x + d, 'linear', NaN);
end

StaticWaveformFamily = struct();
StaticWaveformFamily.method = 'B2_anchored_static_family_plus_diagnostic_peak_alignment';
StaticWaveformFamily.description = ['Peak-aligned xi is used only for static-family ' ...
    'and low-speed localization. The original B2-anchor coordinate is retained.'];
StaticWaveformFamily.sourceStep05File = sourceFile;
StaticWaveformFamily.bladeIds = blades;
StaticWaveformFamily.nominalGapMm = gaps;
StaticWaveformFamily.xiAxisMm = xi;
StaticWaveformFamily.waveformsPeakAlignedMv = Yxi;
StaticWaveformFamily.peakXiMmInB2Coordinate = peakXiMm;
StaticWaveformFamily.peakGridMmInB2Coordinate = peakGridMm;
StaticWaveformFamily.peakVoltageMv = peakVoltageMv;
StaticWaveformFamily.xB2OprAxisMm = x;
StaticWaveformFamily.waveformsB2OprAlignedMv = Yx;
StaticWaveformFamily.waveformsB2FixedPeakAlignedMv = Yb2Fixed;
StaticWaveformFamily.b2ReferenceBlade = 2;
StaticWaveformFamily.b2ReferenceDefinition = 'SG quadratic peak of B2 in each static gap state';
StaticWaveformFamily.staticBaselineMethod = R.staticBaselineMethod;

save(fullfile(outDir, 'StaticWaveformFamily_20251222.mat'), 'StaticWaveformFamily', '-v7.3');
writematrix([blades(:), peakXiMm], fullfile(outDir, 'PeakLocationsByBladeGap_20251222.csv'));

fig = figure('Color','w','Position',[80 80 1380 850]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
colors = turbo(nG);
for ib = 1:nB
    nexttile; hold on; grid on; box on;
    for ig = 1:nG, plot(xi, Yxi(:,ib,ig), 'Color', colors(ig,:), 'LineWidth', .8); end
    xline(0,'k-','LineWidth',.7); title(sprintf('Blade %d: peak-aligned family',blades(ib)));
    xlabel('\xi (mm)'); ylabel('Baseline-subtracted voltage (mV)');
end
saveas(fig, fullfile(outDir, '01_StaticWaveformFamily_PeakAligned_20251222.png'));
fprintf('Static waveform family saved under %s\n', outDir);
