function Plot_EOSearchStrategyStudy()
%PLOT_EOSEARCHSTRATEGYSTUDY Rebuild the publication figure from saved CSV.
root=fileparts(mfilename('fullpath'));outDir=fullfile(root,'output','eo_search_strategy_study');
S=readtable(fullfile(outDir,'eo_search_strategy_summary.csv'),'TextType','string');
snrList=[25 15 5];spans=[9 81];keys=["global_budget","balanced_per_pair","coarse_top_m_multistart"];
labels=["A: global-50","B: stratified-2/pair","C: Top40 x 2"];
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 13]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');cc=[.30 .48 .70;.85 .45 .30;.38 .68 .55];barHandles=[];
for isp=1:numel(spans)
    Q=S(S.smooth_span==spans(isp),:);
    ax=nexttile;Y=metric_matrix(Q,snrList,keys,'eo_correct_rate');b=bar(1:3,100*Y,'grouped');
    for k=1:numel(b),b(k).FaceColor=cc(k,:);end
    set(ax,'XTick',1:3,'XTickLabel',string(snrList));ylim([0 105]);ylabel('Correct EO pair (%)');xlabel('Common SNR (dB)');title(sprintf('SG span = %d',spans(isp)));
    text(ax,-.13,1.04,char('a'+(isp-1)*3),'Units','normalized','FontWeight','bold','FontSize',9);if isp==1,barHandles=b;end
    ax=nexttile;hold on;
    for k=1:numel(keys)
        q=Q.strategy_key==keys(k);[~,o]=ismember(snrList,Q.snr_db(q));v=Q.max_A_error_mm_median(q);v=v(o);
        plot((1:3)+.045*(k-2),1000*v,'-o','Color',cc(k,:),'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor','w');
    end
    set(ax,'XTick',1:3,'XTickLabel',string(snrList));xlabel('Common SNR (dB)');ylabel('Median max amplitude error (\mum)');
    text(ax,-.13,1.04,char('b'+(isp-1)*3),'Units','normalized','FontWeight','bold','FontSize',9);
    ax=nexttile;Y=metric_matrix(Q,snrList,keys,'runtime_s_median');b=bar(1:3,Y,'grouped');
    for k=1:numel(b),b(k).FaceColor=cc(k,:);end
    set(ax,'XTick',1:3,'XTickLabel',string(snrList));ylabel('Median solver time (s)');xlabel('Common SNR (dB)');
    text(ax,-.13,1.04,char('c'+(isp-1)*3),'Units','normalized','FontWeight','bold','FontSize',9);
end
set(findall(fig,'Type','axes'),'TickDir','in','Box','on','FontName','Times New Roman','FontSize',8,'LineWidth',.7);
lg=legend(barHandles,cellstr(labels),'Box','off','Orientation','horizontal');lg.Layout.Tile='south';
sgtitle('EO search allocation controls robustness independently of template smoothing','FontName','Times New Roman','FontSize',10,'FontWeight','bold');
exportgraphics(fig,fullfile(outDir,'eo_search_strategy_comparison.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'eo_search_strategy_comparison.pdf'),'ContentType','vector');
print(fig,fullfile(outDir,'eo_search_strategy_comparison.svg'),'-dsvg');close(fig);
end

function Y=metric_matrix(Q,snrList,keys,field)
Y=nan(numel(snrList),numel(keys));
for i=1:numel(snrList),for k=1:numel(keys),q=Q.snr_db==snrList(i)&Q.strategy_key==keys(k);Y(i,k)=Q.(field)(q);end,end
end
