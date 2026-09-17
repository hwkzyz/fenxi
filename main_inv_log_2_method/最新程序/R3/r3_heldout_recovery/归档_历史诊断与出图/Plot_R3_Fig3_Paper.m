function Files = Plot_R3_Fig3_Paper()
%PLOT_R3_FIG3_PAPER Main-text figure: fixed failure and adaptive rescue.
root=fileparts(mfilename('fullpath')); file=fullfile(root,'output','18_paper_data_sync','r3_paper_result.mat');
if ~isfile(file),error('r3:MissingPaperData','Run Build_R3_PaperDataset first.');end
S=load(file,'Result');R=S.Result;P=R.P;g=sort(P.gMainMm);A=P.amplitudeMainMm;
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 11.8]);
t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t,1); plot_heatmap(ax,R.Summary,"fixed",g,A,'P_vib','Fixed reference'); panel_tag(ax,'a');
ax=nexttile(t,2); plot_heatmap(ax,R.Summary,"adaptive",g,A,'P_vib','Clearance adaptive'); panel_tag(ax,'b');
ax=nexttile(t,3); plot_heatmap(ax,R.Paired,"paired",g,A,'P_rescue','Paired adaptive rescue'); panel_tag(ax,'c');
ax=nexttile(t,4); plot_failure_composition(ax,R.FailureByGap); panel_tag(ax,'d');
colormap(fig,probability_map(256)); cb=colorbar(ax);cb.Layout.Tile='east';cb.Label.String='Probability / fraction';cb.FontName='Times New Roman';cb.FontSize=8;
title(t,'Synchronous-vibration recovery under clearance mismatch','FontName','Times New Roman','FontSize',9,'FontWeight','normal');
outDir=fullfile(root,'output','figures_final');if ~exist(outDir,'dir'),mkdir(outDir);end;base=fullfile(outDir,'Fig3_R3_recovery_paper');
exportgraphics(fig,[base '.png'],'Resolution',600);exportgraphics(fig,[base '.pdf'],'ContentType','vector');exportgraphics(fig,[base '.emf'],'ContentType','vector');
Files=struct('png',[base '.png'],'pdf',[base '.pdf'],'emf',[base '.emf']);
end
function plot_heatmap(ax,T,method,g,A,varName,titleText)
Z=NaN(numel(A),numel(g));for ia=1:numel(A),for ig=1:numel(g),if method=="paired",q=abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;else,q=T.method==method&abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;end;if any(q),Z(ia,ig)=T.(varName)(find(q,1));end;end;end
imagesc(ax,g-P.gReferenceMm,A,Z);set(ax,'YDir','normal','CLim',[0 1]);hold(ax,'on');
for x=(g-P.gReferenceMm(1:end-1)+g(2:end)-P.gReferenceMm(1:end-1))/2 %#ok<NASGU>
end
set(ax,'XTick',g-P.gReferenceMm,'YTick',A);xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'Vibration amplitude, A (mm)');title(ax,titleText);style_axes(ax);
end
function plot_failure_composition(ax,T)
cats=["EO","EO + amplitude","EO + phase","EO + amplitude + phase"];g=sort(unique(T.g_truth_mm));M=zeros(numel(g),numel(cats));
for i=1:numel(g),for j=1:numel(cats),q=T.g_truth_mm==g(i)&T.failure_pattern==cats(j);if any(q),M(i,j)=T.fraction_of_failures(q);end;end;end
bar(ax,g,M,'stacked','BarWidth',.72);ylim(ax,[0 1]);xlabel(ax,'True clearance, g (mm)');ylabel(ax,'Fraction of Fixed failures');title(ax,'Fixed failure-mode composition');legend(ax,cats,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',7);style_axes(ax);
end
function style_axes(ax),set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.75,'TickDir','in','Box','on','Layer','top','XGrid','off','YGrid','off');end
function panel_tag(ax,s),text(ax,-.13,1.06,['(' s ')'],'Units','normalized','FontName','Times New Roman','FontSize',8.5,'HorizontalAlignment','left','VerticalAlignment','bottom');end
function c=probability_map(n),x=[0 .45 .75 1];a=[.16 .23 .32;.32 .61 .70;.94 .76 .33;.72 .18 .20];c=interp1(x,a,linspace(0,1,n),'pchip');c=max(0,min(1,c));end
