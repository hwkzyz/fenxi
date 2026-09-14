function report = Plot_C0_MethodComparisonVerified(outputDir)
%PLOT_C0_METHODCOMPARISONVERIFIED Plot verified Foundation/V1/R5 summaries.
% Values are copied from the frozen comparison summaries recorded in
% R5验证记录_20260914.md. This figure is descriptive only; it does not tune
% EO, frequency, amplitude, or any identification setting.
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', ...
        'diagnostics_20260914', 'C0_method_comparison');
end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

caseName = {'20250527 B1/R01'; '20251222 B1/R01'; '20251222 B5/R04'};
method = {'Foundation'; 'V1 gap-aware'; 'R5'};
% Columns: frequency mean, frequency SD, amplitude mean, amplitude SD.
% 20250527 V1 summary was not part of the accepted formal bundle; retain NaN
% rather than silently reconstructing it.
M = [ ...
    NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; 580.1788,0.0417,NaN,NaN; ...
    582.5353,19.1814,0.1129,0.0054; 582.6485,0.4882,0.1127,0.0038; 582.6463,0.4881,0.1129,0.0038; ...
    633.0219,188.0372,0.1347,0.0588; 631.5735,32.6168,0.1068,0.0232; 631.5717,1.7816,0.1064,0.0219];
caseIdx = [1;1;1;2;2;2;3;3;3];
methodIdx = [1;2;3;1;2;3;1;2;3];
T = table(caseName(caseIdx),method(methodIdx),M(:,1),M(:,2),M(:,3),M(:,4), ...
    'VariableNames',{'case_name','method','frequency_mean_hz','frequency_sd_hz','amplitude_mean_mm','amplitude_sd_mm'});
writetable(T,fullfile(outputDir,'C0_verified_method_summary.csv'));

fig = figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 11]);
t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
colors = [0.25 0.25 0.25; 0.10 0.35 0.75; 0.80 0.15 0.10];
for ip = 1:2
    ax = nexttile(t); hold(ax,'on'); box(ax,'on');
    vals = nan(3,3); errs = nan(3,3);
    for ic=1:3
        q = caseIdx==ic;
        vals(ic,:) = M(q,2*ip-1).'; errs(ic,:) = M(q,2*ip).';
    end
    x = 1:3; dx = [-0.22 0 0.22];
    for im=1:3
        errorbar(ax,x+dx(im),vals(:,im),errs(:,im),'o','Color',colors(im,:), ...
            'MarkerFaceColor',colors(im,:),'LineWidth',1.0,'CapSize',4,'DisplayName',method{im});
    end
    if ip==1, ylabel(ax,'Frequency (Hz)'); else, ylabel(ax,'Amplitude (mm)'); end
    set(ax,'XTick',x,'XTickLabel',{'20250527 B1','20251222 B1','20251222 B5'}, ...
        'FontName','Times New Roman','FontSize',8,'TickDir','in','LineWidth',0.75,'XGrid','off','YGrid','off');
    if ip==1, title(ax,'Verified method comparison (mean +/- SD)','FontWeight','normal'); end
    if ip==2, xlabel(ax,'Test condition'); end
    legend(ax,'Location','best','FontSize',7);
end
exportgraphics(fig,fullfile(outputDir,'C0_verified_method_comparison.png'),'Resolution',300);
close(fig);
report = struct('schema','C0_VERIFIED_METHOD_COMPARISON_V1','outputDir',outputDir, ...
    'status','generated','note','NaN denotes no accepted formal V1 summary for 20250527.');
save(fullfile(outputDir,'C0_VerifiedMethodComparison.mat'),'report','T','-v7.3');
end
