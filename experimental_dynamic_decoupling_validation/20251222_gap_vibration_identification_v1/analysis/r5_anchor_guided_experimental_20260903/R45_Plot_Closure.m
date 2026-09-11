function R45_Plot_Closure()
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results','unified_closure_20260905');
T=readtable(fullfile(out,'UnifiedTransfer_PerAnchor.csv'),'TextType','string');
f=figure('Visible','off','Units','centimeters','Position',[2 2 17 7.5],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for k=1:2
 if k==1, name="R4"; else, name="R5_static"; end
 nexttile; A=T(T.dataset==name & T.method=="PC1" & ~T.hidden_gap,:);
 B=T(T.dataset==name & T.method=="PC1" & T.hidden_gap,:);
 assert(isequal(A(:,{'target_index','anchor_index'}),B(:,{'target_index','anchor_index'})));
 scatter(A.held_surface_nrmse_pct,B.held_surface_nrmse_pct,18,[0 0.45 0.7],'filled');
 hold on; v=[A.held_surface_nrmse_pct;B.held_surface_nrmse_pct]; lim=[max(min(v)*0.7,1e-4),max(v)*1.3];
 plot(lim,lim,'k--','LineWidth',0.7); xlim(lim); ylim(lim); set(gca,'XScale','log','YScale','log');
 xlabel('Known-gap surface NRMSE (%)'); ylabel('Hidden-gap surface NRMSE (%)');
 title(sprintf('(%c) %s',96+k,strrep(char(name),'_',' ')),'FontWeight','normal','Interpreter','none');
 style(gca);
end
savefigures(f,out,'UnifiedTransfer_KnownVsHidden'); close(f);
S=readtable(fullfile(out,'IndependentStrain_PerWindow.csv')); S=S(S.band_low_hz==300,:);
f=figure('Visible','off','Units','centimeters','Position',[2 2 17 8],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for k=1:2
 blades=[1 5]; blade=blades(k); nexttile; a=S(S.blade==blade & S.channel==1,:); b=S(S.blade==blade & S.channel==3,:);
 plot(a.window,a.btt_frequency_hz,'o-','Color',[0.15 0.15 0.15],'MarkerSize',3,'LineWidth',0.9); hold on;
 plot(a.window,a.strain_frequency_hz,'-','Color',[0 0.45 0.7],'LineWidth',1.1);
 plot(b.window,b.strain_frequency_hz,'--','Color',[0.85 0.33 0.1],'LineWidth',1.1);
 xlabel('Window'); ylabel('Frequency (Hz)'); xlim([1 18]);
 title(sprintf('(%c) B%d',96+k,blade),'FontWeight','normal');
 legend({'Main10 + historical R5','Strain AI1-01','Strain AI1-03'},'Location','best','Box','off','FontSize',8);
 style(gca);
end
savefigures(f,out,'IndependentStrain_Frequency'); close(f);
if isfile(fullfile(out,'Ridge_PerWindowSpans.csv'))
 R=readtable(fullfile(out,'Ridge_PerWindowSpans.csv'));
 f=figure('Visible','off','Units','centimeters','Position',[2 2 17 8],'Color','w');
 tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
 nexttile; hold on;
 for blade=[1 5]
  a=R(R.blade==blade,:); plot(a.window,1000*a.frequency_span_hz,'o-','MarkerSize',3,'DisplayName',sprintf('B%d',blade));
 end
 xlabel('Window'); ylabel('Frequency span (mHz)'); title('(a) Five admissible surfaces','FontWeight','normal');
 legend('Location','best','Box','off'); style(gca);
 nexttile; hold on;
 for blade=[1 5]
  a=R(R.blade==blade,:); plot(a.window,1000*a.amplitude_span_mm,'o-','MarkerSize',3,'DisplayName',sprintf('B%d',blade));
 end
 xlabel('Window'); ylabel('Amplitude span (\mum)'); title('(b) Five admissible surfaces','FontWeight','normal');
 legend('Location','best','Box','off'); style(gca);
 savefigures(f,out,'Ridge_DynamicSpans'); close(f);
end
end
function style(ax)
set(ax,'FontName','Times New Roman','FontSize',8,'TickDir','in','Box','on','XGrid','off','YGrid','off');
ax.XLabel.FontSize=9; ax.YLabel.FontSize=9; ax.Title.FontSize=9;
end
function savefigures(f,out,name)
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',300);
exportgraphics(f,fullfile(out,[name '.pdf']),'ContentType','vector');
print(f,fullfile(out,[name '.emf']),'-dmeta');
end
