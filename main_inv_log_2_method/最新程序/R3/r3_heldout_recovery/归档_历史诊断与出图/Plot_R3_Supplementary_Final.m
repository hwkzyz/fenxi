function Files = Plot_R3_Supplementary_Final()
%PLOT_R3_SUPPLEMENTARY_FINAL Null limitation and joint recovery, corrected.
root=fileparts(mfilename('fullpath')); file=fullfile(root,'output','18_paper_data_sync','r3_paper_result.mat');
if ~isfile(file),error('r3:MissingPaperData','Run Build_R3_PaperDataset first.');end
S=load(file,'Result'); R=S.Result; P=R.P; T=R.Summary; g=sort(P.gMainMm); d=g-P.gReferenceMm;
methods=["fixed","adaptive","state_matched_corrected"]; labels=["Fixed reference","Clearance adaptive","State matched (corrected)"];
colors=[.20 .20 .20;0 .45 .62;.78 .29 .20]; marks=['o','s','^']; fig=figure('Color','w','Units','centimeters','Position',[2 2 17 7.2]);
t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact'); ax=nexttile(t); hold(ax,'on');
for i=1:3,y=extract_line(T,methods(i),g,0,'P_FP');plot(ax,d,y,'-','Color',colors(i,:),'Marker',marks(i),'LineWidth',1.1,'MarkerSize',4,'DisplayName',labels(i));end
xlabel(ax,'Clearance mismatch, \Delta g (mm)'); ylabel(ax,'False-positive probability, P_{FP}'); ylim(ax,[0 1]); legend(ax,'Location','best','Box','off'); panel_tag(ax,'a'); style_axes(ax);
ax=nexttile(t); hold(ax,'on');
for i=1:3
 y=NaN(size(g)); for j=1:numel(g),q=T.method==methods(i)&T.g_truth_mm==g(j)&T.A_true_mm>0;y(j)=mean(T.P_joint(q),'omitnan');end
 plot(ax,d,y,'-','Color',colors(i,:),'Marker',marks(i),'LineWidth',1.1,'MarkerSize',4,'DisplayName',labels(i));
end
yline(ax,P.targetProbability,'--','0.90 criterion','LabelHorizontalAlignment','left','Color',[.35 .35 .35],'FontName','Times New Roman','FontSize',8);
xlabel(ax,'Clearance mismatch, \Delta g (mm)'); ylabel(ax,'Mean joint-success probability, P_{joint}'); ylim(ax,[0 1]); panel_tag(ax,'b'); style_axes(ax);
outDir=fullfile(root,'output','figures_final'); if ~exist(outDir,'dir'),mkdir(outDir);end; base=fullfile(outDir,'FigS_R3_null_and_joint_corrected');
exportgraphics(fig,[base '.png'],'Resolution',600); exportgraphics(fig,[base '.pdf'],'ContentType','vector'); exportgraphics(fig,[base '.emf'],'ContentType','vector');
Files=struct('png',[base '.png'],'pdf',[base '.pdf'],'emf',[base '.emf']);
end
function y=extract_line(T,method,g,A,varName),y=NaN(size(g));for i=1:numel(g),q=T.method==method&abs(T.g_truth_mm-g(i))<1e-12&abs(T.A_true_mm-A)<1e-12;if any(q),y(i)=T.(varName)(find(q,1));end;end;end
function style_axes(ax),set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.75,'TickDir','in','Box','on','XGrid','off','YGrid','off');end
function panel_tag(ax,s),text(ax,-.13,1.06,['(' s ')'],'Units','normalized','FontName','Times New Roman','FontSize',8.5,'HorizontalAlignment','left','VerticalAlignment','bottom');end
