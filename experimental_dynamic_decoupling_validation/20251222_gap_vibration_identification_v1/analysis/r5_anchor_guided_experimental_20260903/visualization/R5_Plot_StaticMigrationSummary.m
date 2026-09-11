function file = R5_Plot_StaticMigrationSummary(T, mode, file)
%R5_PLOT_STATICMIGRATIONSUMMARY Method-level static migration comparison.
fig=figure('Visible','off','Color','w'); tiledlayout(1,2,'TileSpacing','compact');
nexttile; boxchart(categorical(T.blade_id),T.relative_rmse); grid on;
xlabel('blade'); ylabel('relative RMSE'); title(['R5 ',strrep(mode,'_',' ')]);
nexttile; scatter(T.nominal_gap_mm,T.peak_x_error_mm,24,T.blade_id,'filled'); grid on;
xlabel('nominal static gap (mm)'); ylabel('peak x error (mm)'); colorbar;
exportgraphics(fig,file,'Resolution',180); close(fig);
end
