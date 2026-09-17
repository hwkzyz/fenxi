function Main_05_PC1Figures()
%MAIN_05_PC1FIGURES Plot saved PC1/open-set results only.
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, 'output', 'surface_identification_gap_law');
pc1 = load(fullfile(outDir, 'PredictedSurfaceFamilyLocalizationResults.mat'));
open = load(fullfile(outDir, 'OpenSetSurfaceTransferGapLawResults.mat'));
fig = figure('Color','w','Position',[100 100 1050 700]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; bar(pc1.summary.target_tilt_deg,pc1.summary.mean_transfer_NRMSE_pct);
xlabel('Hidden tilt (deg)'); ylabel('PC1 transfer NRMSE (%)'); grid on; box on;
nexttile; bar(open.summary.target_tilt_deg,open.summary.mean_transfer_NRMSE_pct);
xlabel('Hidden tilt (deg)'); ylabel('Open-set transfer NRMSE (%)'); grid on; box on;
nexttile; plot(pc1.T.anchorGap,pc1.T.transferNrmsePct,'.','MarkerSize',8);
xlabel('Anchor gap (mm)'); ylabel('PC1 transfer NRMSE (%)'); grid on; box on;
nexttile;
if ismember('pc1_explained_ratio',pc1.overall.Properties.VariableNames)
    bar([pc1.overall.pc1_explained_ratio,pc1.overall.extrapolation_fraction]);
    ylim([0 1]); set(gca,'XTickLabel',{'PC1 explained','Extrapolation'}); ylabel('Fraction'); grid on; box on;
else
    axis off; text(0.1,0.5,'Run Main_02 again for diagnostics','FontSize',10);
end
exportgraphics(fig,fullfile(outDir,'Fig_PC1_Validation_Summary.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig_PC1_Validation_Summary.pdf'),'ContentType','vector');
close(fig);
fprintf('Figures written from saved result files to:\n%s\n',outDir);
end
