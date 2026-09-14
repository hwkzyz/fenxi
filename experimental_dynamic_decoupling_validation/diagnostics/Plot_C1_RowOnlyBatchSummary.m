function report = Plot_C1_RowOnlyBatchSummary(resultFile, outputDir)
%PLOT_C1_ROWONLYBATCHSUMMARY Batch audit for the public out.rows schema.
% This is intentionally separate from Plot_C_DynamicBatchSummary, which
% requires the larger Result.WindowResult schema.
if nargin<2 || isempty(outputDir), outputDir=fullfile(fileparts(resultFile),'batch'); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
S=load(resultFile); assert(isfield(S,'out') && isfield(S.out,'rows'), ...
    'C1:RowOnly','Expected out.rows in %s.',resultFile);
r=S.out.rows; n=numel(r); w=(1:n).';
getv=@(name) localVector(r,name,n);
f=getv('frequency_hz'); a=getv('amplitude_mm'); rm=getv('rmse_mv');
sp=getv('support_fraction'); eo=getv('EO');
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 13]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(w,f,'o-','LineWidth',1); ylabel('Frequency (Hz)'); title('R5 frequency','FontWeight','normal');
nexttile; plot(w,a,'o-','LineWidth',1); ylabel('Amplitude (mm)'); title('R5 amplitude','FontWeight','normal');
nexttile; plot(w,rm,'o-','LineWidth',1); ylabel('RMSE (mV)'); title('Full-wave residual','FontWeight','normal');
nexttile; stairs(w,sp,'LineWidth',1); hold on; plot(w,eo,'x','LineWidth',1); ylabel('Support / EO'); title('Support fraction and EO','FontWeight','normal'); legend('support','EO','Location','best');
for k=1:4, ax=nexttile(k); xlabel(ax,'Window'); box(ax,'on'); set(ax,'FontName','Times New Roman','FontSize',8,'TickDir','in','XGrid','off','YGrid','off'); end
exportgraphics(fig,fullfile(outputDir,'C1_rowonly_batch_summary.png'),'Resolution',220); close(fig);
T=table(w,f,a,rm,sp,eo,'VariableNames',{'window_id','frequency_hz','amplitude_mm','rmse_mv','support_fraction','EO'});
writetable(T,fullfile(outputDir,'C1_rowonly_batch_summary.csv'));
report=struct('schema','C1_ROWONLY_BATCH_SUMMARY_V1','nWindows',n,'outputDir',outputDir);
save(fullfile(outputDir,'C1_RowOnlyBatchSummary.mat'),'report','T','-v7.3');
end

function v=localVector(r,name,n)
v=nan(n,1);
for i=1:n
    if isfield(r(i),name)
        x=r(i).(name); if isnumeric(x) && isscalar(x), v(i)=double(x); end
    end
end
end
