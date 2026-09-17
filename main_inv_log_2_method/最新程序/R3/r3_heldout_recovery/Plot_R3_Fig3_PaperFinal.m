function Files=Plot_R3_Fig3_PaperFinal()
root=fileparts(mfilename('fullpath')); S=load(fullfile(root,'output','18_paper_data_sync','r3_paper_result.mat'),'Result'); R=S.Result; P=R.P; g=sort(P.gMainMm); A=P.amplitudeMainMm;
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 11.8]); t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t,1); drawmap(ax,R.Summary,"fixed",g,A,'P_vib','Fixed reference'); tag(ax,'a');
ax=nexttile(t,2); drawmap(ax,R.Summary,"adaptive",g,A,'P_vib','Clearance adaptive'); tag(ax,'b');
ax=nexttile(t,3); drawmap(ax,R.Paired,"paired",g,A,'P_rescue','Paired adaptive rescue'); tag(ax,'c');
ax=nexttile(t,4); drawfail(ax,R.FailureByGap); tag(ax,'d');
colormap(fig,parula(256)); cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String='Probability / fraction'; cb.FontName='Times New Roman'; cb.FontSize=8;
title(t,'Synchronous-vibration recovery under clearance mismatch','FontName','Times New Roman','FontSize',9,'FontWeight','normal');
out=fullfile(root,'output','figures_final'); if ~exist(out,'dir'),mkdir(out);end; base=fullfile(out,'Fig3_R3_recovery_paper_final');
exportgraphics(fig,[base '.png'],'Resolution',600); exportgraphics(fig,[base '.pdf'],'ContentType','vector'); exportgraphics(fig,[base '.emf'],'ContentType','vector'); Files=struct('png',[base '.png'],'pdf',[base '.pdf'],'emf',[base '.emf']);
end
function drawmap(ax,T,m,g,A,v,ttl)
Z=NaN(numel(A),numel(g));
for ia=1:numel(A)
 for ig=1:numel(g)
  if m=="paired", q=abs(T.g_truth_mm-g(ig))<1e-12 & abs(T.A_true_mm-A(ia))<1e-12;
  else, q=T.method==m & abs(T.g_truth_mm-g(ig))<1e-12 & abs(T.A_true_mm-A(ia))<1e-12; end
  if any(q), Z(ia,ig)=T.(v)(find(q,1)); end
 end
end
imagesc(ax,g-0.8,A,Z); set(ax,'YDir','normal','CLim',[0 1],'XTick',g-0.8,'YTick',A); xlabel(ax,'Clearance mismatch, \Delta g (mm)'); ylabel(ax,'Vibration amplitude, A (mm)'); title(ax,ttl); style(ax);
end
function drawfail(ax,T)
cats=["EO","EO + amplitude","EO + phase","EO + amplitude + phase"]; gg=sort(unique(T.g_truth_mm)); M=zeros(numel(gg),numel(cats));
for i=1:numel(gg)
 for j=1:numel(cats)
  q=T.g_truth_mm==gg(i)&T.failure_pattern==cats(j); if any(q), M(i,j)=T.fraction_of_failures(find(q,1)); end
 end
end
bar(ax,gg,M,'stacked'); ylim(ax,[0 1]); xlabel(ax,'True clearance, g (mm)'); ylabel(ax,'Fraction of Fixed failures'); title(ax,'Fixed failure-mode composition'); legend(ax,cats,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',7); style(ax);
end
function style(ax),set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.75,'TickDir','in','Box','on','Layer','top','XGrid','off','YGrid','off');end
function tag(ax,s),text(ax,-.13,1.06,['(' s ')'],'Units','normalized','FontName','Times New Roman','FontSize',8.5);end
