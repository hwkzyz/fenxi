function report = Plot_A4_SensorSensitivityAudit(csvFile,outputDir)
%PLOT_A4_SENSORSENSITIVITYAUDIT Visualize per-sensor static gap sensitivity.
if nargin<1||isempty(csvFile), csvFile=fullfile(fileparts(mfilename('fullpath')),'..','..','fenxi_release_20260911','20250527_gap_vibration_identification_v1','latest_programs_20260907','results','r5_sensor_sensitivity_audit.csv'); end
if nargin<2||isempty(outputDir), outputDir=fullfile(fileparts(mfilename('fullpath')),'..','results','diagnostics_20260914','A4_sensitivity'); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
T=readtable(csvFile); fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 10]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile; bar(ax,T.sensor_id,[T.median_abs_dFdg T.p90_abs_dFdg]); box(ax,'on'); ylabel(ax,'|F_g| (mV/mm)'); xlabel(ax,'Sensor'); legend(ax,'median','p90','Location','best');
xticks(ax,T.sensor_id); xticklabels(ax,compose('S%d',T.sensor_id));
ax=nexttile; bar(ax,T.sensor_id,T.anchor_rmse_mv); box(ax,'on'); ylabel(ax,'Anchor replay RMSE (mV)'); xlabel(ax,'Sensor');
xticks(ax,T.sensor_id); xticklabels(ax,compose('S%d',T.sensor_id));
sgtitle(fig,'20250527 | A4 Static Sensitivity Audit','FontName','Times New Roman','FontSize',9,'FontWeight','normal');
for ax=findall(fig,'Type','axes')', set(ax,'FontName','Times New Roman','FontSize',8,'TickDir','in','LineWidth',.75,'XGrid','off','YGrid','off'); end
exportgraphics(fig,fullfile(outputDir,'A4_static_sensitivity_audit.png'),'Resolution',300);
exportgraphics(fig,fullfile(outputDir,'A4_static_sensitivity_audit.pdf'),'ContentType','vector'); close(fig);
report=struct('schema','A4_STATIC_SENSITIVITY_AUDIT_V1','source',csvFile,'status','generated');
save(fullfile(outputDir,'A4_StaticSensitivityAudit.mat'),'report','T','-v7.3');
end
