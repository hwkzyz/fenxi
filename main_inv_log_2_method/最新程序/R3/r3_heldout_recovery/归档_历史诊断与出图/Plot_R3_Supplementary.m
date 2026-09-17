function Files = Plot_R3_Supplementary(resultFile)
%PLOT_R3_SUPPLEMENTARY Plot null false positives and joint recovery evidence.

root=fileparts(mfilename('fullpath'));
if nargin<1||isempty(resultFile)
    resultFile=fullfile(root,'output','05_main_15db','main_result.mat');
end
if ~exist(resultFile,'file')
    error('r3:MissingMainResult','Run Run_05_Main15dB before plotting supplementary results.');
end
S=load(resultFile,'Result');R=S.Result;P=R.P;T=R.Summary;
g=sort(P.gMainMm);d=g-P.gReferenceMm;methods=["fixed","adaptive","state_matched"];
labels=["Fixed reference","Clearance adaptive","State matched"];
colors=[0.20 0.20 0.20;0.00 0.45 0.62;0.78 0.29 0.20];marks=['o','s','^'];

fig=figure('Color','w','Units','centimeters','Position',[2 2 17 7.2]);
t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t);hold(ax,'on');
for i=1:3
    y=extract_line(T,methods(i),g,0,'P_FP');
    plot(ax,d,y,'-','Color',colors(i,:),'Marker',marks(i),'LineWidth',1.1,...
        'MarkerSize',4,'DisplayName',labels(i));
end
xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'False-positive probability, P_{FP}');
ylim(ax,[0 1]);legend(ax,'Location','best','Box','off');panel_tag(ax,'a');style_axes(ax);

ax=nexttile(t);hold(ax,'on');
for i=1:3
    q=T.method==methods(i)&T.A_true_mm>0;
    G=groupsummary(T(q,:),{'g_truth_mm'},'mean','P_joint');
    [gx,k]=sort(G.g_truth_mm);plot(ax,gx-P.gReferenceMm,G.mean_P_joint(k),'-',...
        'Color',colors(i,:),'Marker',marks(i),'LineWidth',1.1,'MarkerSize',4,...
        'DisplayName',labels(i));
end
yline(ax,P.targetProbability,'--','0.90 criterion','LabelHorizontalAlignment','left',...
    'Color',[0.35 0.35 0.35],'FontName','Times New Roman','FontSize',8);
xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'Mean joint-success probability, P_{joint}');
ylim(ax,[0 1]);panel_tag(ax,'b');style_axes(ax);

outDir=fullfile(root,'output','figures');if ~exist(outDir,'dir'),mkdir(outDir);end
base=fullfile(outDir,'FigS_R3_false_positive_and_joint');
exportgraphics(fig,[base '.png'],'Resolution',600);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
exportgraphics(fig,[base '.emf'],'ContentType','vector');
Files=struct('png',[base '.png'],'pdf',[base '.pdf'],'emf',[base '.emf']);
end

function y=extract_line(T,method,g,A,varName)
y=NaN(size(g));
for i=1:numel(g)
    q=T.method==method&abs(T.g_truth_mm-g(i))<1e-12&abs(T.A_true_mm-A)<1e-12;
    if any(q),y(i)=T.(varName)(find(q,1));end
end
end
function style_axes(ax)
set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',0.75,...
    'TickDir','in','Box','on','XGrid','off','YGrid','off');
end
function panel_tag(ax,s)
text(ax,-0.13,1.06,['(' s ')'],'Units','normalized','FontName','Times New Roman',...
    'FontSize',8.5,'HorizontalAlignment','left','VerticalAlignment','bottom');
end
