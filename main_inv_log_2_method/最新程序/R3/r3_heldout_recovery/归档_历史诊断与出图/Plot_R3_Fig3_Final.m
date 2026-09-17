function Files = Plot_R3_Fig3_Final()
%PLOT_R3_FIG3_FINAL Paper-facing R3 figure with corrected state-matched data.
root=fileparts(mfilename('fullpath')); file=fullfile(root,'output','18_paper_data_sync','r3_paper_result.mat');
if ~isfile(file),error('r3:MissingPaperData','Run Build_R3_PaperDataset first.');end
S=load(file,'Result'); R=S.Result; P=R.P; g=sort(P.gMainMm); A=P.amplitudeMainMm;
methods=["fixed","adaptive","state_matched_corrected"];
titles=["Fixed reference","Clearance adaptive","State matched (corrected)"];
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 11.8]);
t=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for i=1:3
 ax=nexttile(t,i); Z=domain_matrix(R.Summary,methods(i),g,A,'P_vib'); imagesc(ax,g-P.gReferenceMm,A,Z); set(ax,'YDir','normal','CLim',[0 1]);
 hold(ax,'on'); finiteZ=Z(isfinite(Z)); if ~isempty(finiteZ)&&min(finiteZ)<=P.targetProbability&&max(finiteZ)>=P.targetProbability,contour(ax,g-P.gReferenceMm,A,Z,[P.targetProbability P.targetProbability],'LineColor',[.12 .12 .12],'LineWidth',1.1,'LineStyle','--');end
 title(ax,titles(i)); xlabel(ax,'Clearance mismatch, \Delta g (mm)'); ylabel(ax,'Vibration amplitude, A (mm)'); panel_tag(ax,char('a'+i-1)); style_axes(ax);
end
ax=nexttile(t,4); Z=paired_matrix(R.Paired,g,A,'P_rescue'); imagesc(ax,g-P.gReferenceMm,A,Z); set(ax,'YDir','normal','CLim',[0 1]);
title(ax,'Adaptive rescue'); xlabel(ax,'Clearance mismatch, \Delta g (mm)'); ylabel(ax,'Vibration amplitude, A (mm)'); panel_tag(ax,'d'); style_axes(ax);
colormap(fig,probability_map(256)); cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String='Probability'; cb.FontName='Times New Roman'; cb.FontSize=8;
title(t,'Synchronous-vibration recovery under held-out clearance states','FontName','Times New Roman','FontSize',9,'FontWeight','normal');
outDir=fullfile(root,'output','figures_final'); if ~exist(outDir,'dir'),mkdir(outDir);end; base=fullfile(outDir,'Fig3_R3_recovery_corrected');
exportgraphics(fig,[base '.png'],'Resolution',600); exportgraphics(fig,[base '.pdf'],'ContentType','vector'); exportgraphics(fig,[base '.emf'],'ContentType','vector');
Files=struct('png',[base '.png'],'pdf',[base '.pdf'],'emf',[base '.emf']);
end
function Z=domain_matrix(T,method,g,A,varName)
Z=NaN(numel(A),numel(g)); for ia=1:numel(A),for ig=1:numel(g),q=T.method==method&abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;if any(q),Z(ia,ig)=T.(varName)(find(q,1));end;end;end
end
function Z=paired_matrix(T,g,A,varName)
Z=NaN(numel(A),numel(g)); for ia=1:numel(A),for ig=1:numel(g),q=abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;if any(q),Z(ia,ig)=T.(varName)(find(q,1));end;end;end
end
function style_axes(ax),set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.75,'TickDir','in','Box','on','Layer','top','XGrid','off','YGrid','off');end
function panel_tag(ax,s),text(ax,-.13,1.06,['(' s ')'],'Units','normalized','FontName','Times New Roman','FontSize',8.5,'HorizontalAlignment','left','VerticalAlignment','bottom');end
function c=probability_map(n),x=[0 .45 .75 1];a=[.16 .23 .32;.32 .61 .70;.94 .76 .33;.72 .18 .20];c=interp1(x,a,linspace(0,1,n),'pchip');c=max(0,min(1,c));end
