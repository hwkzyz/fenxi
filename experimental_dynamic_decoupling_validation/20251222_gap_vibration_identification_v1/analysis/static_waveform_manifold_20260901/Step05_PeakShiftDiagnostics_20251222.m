%% 05: Diagnose peak position versus nominal gap in the B2 coordinate
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath')); outDir = fullfile(thisDir,'results');
S = load(fullfile(outDir,'StaticWaveformFamily_20251222.mat'),'StaticWaveformFamily'); F=S.StaticWaveformFamily;
g = F.nominalGapMm(:); p = F.peakXiMmInB2Coordinate;
nB = numel(F.bladeIds); rows = []; summary = nan(nB,5);
for ib=1:nB
    ok=isfinite(p(ib,:)).';
    if nnz(ok)>=3
        c=polyfit(g(ok),p(ib,ok).',1); pred=polyval(c,g(ok));
        slope=c(1); rmse=sqrt(mean((pred-p(ib,ok).').^2));
    else
        slope=NaN; rmse=NaN;
    end
    summary(ib,:)=[F.bladeIds(ib),slope,range(p(ib,ok)),std(p(ib,ok)),rmse];
    rows=[rows; [repmat(F.bladeIds(ib),nnz(ok),1),g(ok),p(ib,ok).']]; %#ok<AGROW>
end
writematrix(rows,fullfile(outDir,'05_PeakPositionByGap_20251222.csv'));
writematrix(summary,fullfile(outDir,'05_PeakShiftSummary_20251222.csv'));
fig=figure('Color','w','Position',[120 120 1000 620]); hold on; grid on; box on;
cc=lines(nB);
for ib=1:nB, plot(g,p(ib,:),'-o','Color',cc(ib,:),'LineWidth',1.1,'MarkerSize',4,'DisplayName',sprintf('B%d',F.bladeIds(ib))); end
xlabel('B2 nominal gap (mm)'); ylabel('SG/quadratic peak x in B2 coordinate (mm)');
title('Peak position retained as a static-state diagnostic'); legend('Location','best');
saveas(fig,fullfile(outDir,'05_PeakPositionVsGap_20251222.png'));
fprintf('Peak-shift diagnostics saved. B2 range = %.6g mm.\n',summary(F.bladeIds==2,3));

