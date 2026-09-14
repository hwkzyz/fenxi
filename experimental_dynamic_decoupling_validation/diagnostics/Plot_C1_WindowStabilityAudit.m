function report = Plot_C1_WindowStabilityAudit(csvFiles, labels, outputDir)
%PLOT_C1_WINDOWSTABILITYAUDIT Plot window-level stability from comparison CSVs.
% csvFiles are R5_FormalMethodComparison_*_windows.csv files. The plot is
% read-only and exposes frequency, amplitude, EO agreement, and RMSE change.
if nargin < 1 || isempty(csvFiles)
    root = fullfile(fileparts(mfilename('fullpath')),'..','..','fenxi_release_20260911', ...
        '20251222_gap_vibration_identification_v1','latest_programs_20260905','results','gap_aware','comparison');
    csvFiles = {fullfile(root,'R5_FormalMethodComparison_B1_S123_R01_windows.csv'), ...
                fullfile(root,'R5_FormalMethodComparison_B5_S123_R04_windows.csv')};
end
if nargin < 2 || isempty(labels), labels = {'20251222 B1/R01','20251222 B5/R04'}; end
if nargin < 3 || isempty(outputDir), outputDir = fullfile(fileparts(mfilename('fullpath')),'..','results','diagnostics_20260914','C1_window_audit'); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
fig = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 14]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
rows = table();
for k=1:numel(csvFiles)
    T = readtable(csvFiles{k});
    if isempty(rows), rows=T; rows.case_name=repmat(string(labels{k}),height(T),1); else, T.case_name=repmat(string(labels{k}),height(T),1); rows=[rows;T]; end %#ok<AGROW>
    x=T.window_id;
    ax=nexttile(1); hold(ax,'on'); plot(ax,x,T.proposed_frequency_hz,'o-','DisplayName',labels{k});
    ax=nexttile(2); hold(ax,'on'); plot(ax,x,T.proposed_amplitude_mm,'o-','DisplayName',labels{k});
    ax=nexttile(3); hold(ax,'on'); plot(ax,x,T.rmse_reduction_percent,'o-','DisplayName',labels{k});
    ax=nexttile(4); hold(ax,'on'); stairs(ax,x,T.eo_agreement,'LineWidth',1.1,'DisplayName',labels{k});
end
titles={'R5 frequency by window','R5 amplitude by window','RMSE reduction vs Foundation','EO agreement (1 = agree)'};
ylabels={'Frequency (Hz)','Amplitude (mm)','RMSE reduction (%)','EO agreement'};
for i=1:4
    ax=nexttile(i); box(ax,'on'); set(ax,'FontName','Times New Roman','FontSize',8,'TickDir','in','LineWidth',.75,'XGrid','off','YGrid','off');
    title(ax,titles{i},'FontWeight','normal'); ylabel(ax,ylabels{i}); xlabel(ax,'Window'); legend(ax,'Location','best','FontSize',7);
end
exportgraphics(fig,fullfile(outputDir,'C1_window_stability_audit.png'),'Resolution',300); close(fig);
writetable(rows,fullfile(outputDir,'C1_window_stability_audit.csv'));
report=struct('schema','C1_WINDOW_STABILITY_AUDIT_V1','outputDir',outputDir,'status','generated','sourceFiles',{csvFiles});
save(fullfile(outputDir,'C1_WindowStabilityAudit.mat'),'report','rows','-v7.3');
end
