function C = R5_Compare_Main10_Waveforms(baselineFile,r5File,outputDir,tag)
% Compare frozen Main10 and R5-sidecar results without changing either run.
if nargin<4 || isempty(tag), tag='comparison'; end
if ~isfolder(outputDir), mkdir(outputDir); end
A=load(baselineFile,'Result','Trend','Summary');
B=load(r5File,'Result','Trend','Summary');
assert(numel(A.Result.WindowResult)==numel(B.Result.WindowResult), ...
    'R5:WindowCountMismatch','Baseline and R5 window counts differ.');
n=numel(A.Result.WindowResult);
window_id=(1:n).'; baseline_eo=nan(n,1); r5_eo=nan(n,1);
baseline_frequency_hz=nan(n,1); r5_frequency_hz=nan(n,1);
baseline_rmse_mv=nan(n,1); r5_rmse_mv=nan(n,1);
for k=1:n
    fa=A.Result.WindowResult(k).modelFits.gap_only;
    fb=B.Result.WindowResult(k).modelFits.gap_only;
    baseline_eo(k)=fa.EO; r5_eo(k)=fb.EO;
    baseline_frequency_hz(k)=fa.freqHz; r5_frequency_hz(k)=fb.freqHz;
    baseline_rmse_mv(k)=fa.plainRmseMv; r5_rmse_mv(k)=fb.plainRmseMv;
end
C=table(window_id,baseline_eo,r5_eo,baseline_frequency_hz,r5_frequency_hz, ...
    r5_frequency_hz-baseline_frequency_hz,baseline_rmse_mv,r5_rmse_mv, ...
    r5_rmse_mv-baseline_rmse_mv, ...
    'VariableNames',{'window_id','baseline_eo','r5_eo','baseline_frequency_hz', ...
    'r5_frequency_hz','frequency_difference_hz','baseline_rmse_mv','r5_rmse_mv', ...
    'rmse_difference_mv'});
writetable(C,fullfile(outputDir,sprintf('R5_Main10_WindowComparison_%s.csv',tag)));

% Use the R5 best-fit window and compare both predictions on the identical
% measured point bundle.  This is a direct waveform check, not just a
% frequency-summary comparison.
k=B.Result.BestWindowIndex;
wa=A.Result.WindowResult(k); wb=B.Result.WindowResult(k);
assert(isequaln(wa.bundle.sensorIndex,wb.bundle.sensorIndex) && ...
    max(abs(wa.bundle.X-wb.bundle.X),[],'all')<1e-12 && ...
    max(abs(wa.bundle.V-wb.bundle.V),[],'all')<1e-9, ...
    'R5:PointBundleMismatch','Baseline and R5 best-window point bundles differ.');
sensorIds=wb.bundle.sensorIds(:).';
f=figure('Visible','off','Color','w','Position',[100 100 1200 320*numel(sensorIds)]);
tiledlayout(numel(sensorIds),1,'TileSpacing','compact','Padding','compact');
for j=1:numel(sensorIds)
    nexttile; idx=wb.bundle.sensorIndex==j;
    [xs,ord]=sort(wb.bundle.X(idx)); vv=wb.bundle.V(idx);
    pa=wa.modelFits.gap_only.VPred(idx); pb=wb.modelFits.gap_only.VPred(idx);
    plot(xs,vv(ord),'k.','MarkerSize',4,'DisplayName','measured'); hold on; grid on; box on;
    plot(xs,pa(ord),'Color',[0.25 0.45 0.85],'LineWidth',1.0,'DisplayName','baseline Main10');
    plot(xs,pb(ord),'Color',[0.85 0.25 0.20],'LineWidth',1.0,'DisplayName','R5 sidecar');
    ylabel(sprintf('CH%d / mV',sensorIds(j)));
    if j==1
        title(sprintf('%s, window %d: measured and reconstructed waveforms',tag,k),'Interpreter','none');
        legend('Location','best');
    end
end
xlabel('x / mm');
exportgraphics(f,fullfile(outputDir,sprintf('R5_Main10_WaveformComparison_%s.png',tag)),'Resolution',200);
close(f);

f=figure('Visible','off','Color','w');
plot(window_id,baseline_frequency_hz,'o-','LineWidth',1.1,'DisplayName','baseline Main10'); hold on; grid on; box on;
plot(window_id,r5_frequency_hz,'s-','LineWidth',1.1,'DisplayName','R5 sidecar');
xlabel('window'); ylabel('frequency / Hz'); title(sprintf('%s frequency comparison',tag),'Interpreter','none'); legend('Location','best');
exportgraphics(f,fullfile(outputDir,sprintf('R5_Main10_FrequencyComparison_%s.png',tag)),'Resolution',200); close(f);
end
