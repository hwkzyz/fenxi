function Files=Plot_R3_Fig3_SixPanel()
%PLOT_R3_FIG3_SIXPANEL Six-panel evidence-chain figure for the main text.
root=fileparts(mfilename('fullpath')); dataFile=fullfile(root,'output','18_paper_data_sync','r3_paper_result.mat');
if ~isfile(dataFile),error('r3:MissingPaperData','Run Build_R3_PaperDataset first.');end
S=load(dataFile,'Result'); R=S.Result; P=R.P; g=sort(P.gMainMm); A=P.amplitudeMainMm;
fig=figure('Color','w','Units','centimeters','Position',[1 1 17 17.2]); t=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t,1); draw_protocol(ax,P); tag(ax,'a');
ax=nexttile(t,2); draw_representative(ax,root,P); tag(ax,'b');
ax=nexttile(t,3); draw_map(ax,R.Summary,"fixed",g,A,'P_vib','Fixed reference'); tag(ax,'c');
ax=nexttile(t,4); draw_map(ax,R.Summary,"adaptive",g,A,'P_vib','Clearance adaptive'); tag(ax,'d');
ax=nexttile(t,5); draw_map(ax,R.Paired,"paired",g,A,'P_rescue','Paired adaptive rescue'); tag(ax,'e');
ax=nexttile(t,6); draw_fail(ax,R.FailureByGap); tag(ax,'f');
colormap(fig,parula(256)); title(t,'Synchronous-vibration recovery under matched and held-out clearance states','FontName','Times New Roman','FontSize',9,'FontWeight','normal');
out=fullfile(root,'output','figures_final');if ~exist(out,'dir'),mkdir(out);end;base=fullfile(out,'Fig3_R3_six_panel_evidence');
exportgraphics(fig,[base '.png'],'Resolution',600);exportgraphics(fig,[base '.pdf'],'ContentType','vector');exportgraphics(fig,[base '.emf'],'ContentType','vector');Files=struct('png',[base '.png'],'pdf',[base '.pdf'],'emf',[base '.emf']);
end

function draw_protocol(ax,P)
hold(ax,'on'); allg=P.gAllMm; plot(ax,allg,zeros(size(allg)),'-','Color',[.75 .75 .75],'LineWidth',.8);
cal=P.gCalMm; op=P.gOperatingMm; scatter(ax,cal,zeros(size(cal)),34,[.25 .25 .25],'filled'); scatter(ax,op,zeros(size(op)),38,[0 .45 .62],'filled'); scatter(ax,P.gReferenceMm,zeros(size(P.gReferenceMm)),64,[.78 .29 .20],'o','LineWidth',1.2);
text(ax,P.gReferenceMm,.11,'g_{ref}','HorizontalAlignment','center','FontSize',8);text(ax,.33,.23,'calibration states','Color',[.2 .2 .2],'FontSize',8);text(ax,.98,-.20,'held-out operating states','Color',[0 .35 .5],'FontSize',8,'HorizontalAlignment','center');
xlim(ax,[.15 1.55]);ylim(ax,[-.32 .32]);set(ax,'YTick',[],'XTick',.2:.1:1.5);xlabel(ax,'Clearance, g (mm)');title(ax,'Frozen calibration and held-out test design');style(ax);
text(ax,.22,-.28,'raw high-speed curve \rightarrow noisy record \rightarrow calibration-only inverse model','FontSize',7.4,'Interpreter','tex');
end

function draw_representative(ax,root,P)
mainDir=fileparts(root);addpath(fullfile(root,'local_func'),'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
C=load(fullfile(root,'output','02_reference_calibration','r3_reference_calibration_C0.mat'),'Result'); C0=C.Result.C0; sigmaV=C.Result.sigmaV; R3=r3_build_context(mainDir,P);
g=.7; A=.20; rep=1; seed=P.mainHighSeeds(rep); [~,map,truth]=r3_generate_high_case(R3,g,A,P.phaseValuesRad(rep),sigmaV,seed); F=r3_run_three_methods_global_final(R3,C0,map,g);
t=map.t_v(:); y=map.V_a(:); n=min([numel(t),numel(y),numel(F.fixed.fit.VFit),numel(F.adaptive.fit.VFit),2500]);
plot(ax,t(1:n),y(1:n),'.','Color',[.35 .35 .35],'MarkerSize',2,'DisplayName','Measured');hold(ax,'on');plot(ax,t(1:n),F.fixed.fit.VFit(1:n),'-','Color',[.2 .2 .2],'LineWidth',.8,'DisplayName','Fixed fit');plot(ax,t(1:n),F.adaptive.fit.VFit(1:n),'-','Color',[0 .45 .62],'LineWidth,.9,'DisplayName','Adaptive fit');
xlabel(ax,'Time (s)');ylabel(ax,'Sensor voltage (V)');title(ax,'Representative paired record: g=0.7 mm, A=0.20 mm');legend(ax,'Location','best','Box','off','FontSize',7);style(ax);
txt=sprintf('truth: EO=10, A=%.2f mm\nFixed: EO=%d, A=%.3f mm\nAdaptive: EO=%d, A=%.3f mm, g=%.3f mm',A,F.fixed.fit.eo_id,F.fixed.fit.A_id,F.adaptive.fit.eo_id,F.adaptive.fit.A_id,F.adaptive.fit.g_used);text(ax,.02,.97,txt,'Units','normalized','VerticalAlignment','top','FontSize',7.2,'BackgroundColor','w','Margin',2);
end

function draw_map(ax,T,m,g,A,v,ttl)
Z=NaN(numel(A),numel(g));for ia=1:numel(A),for ig=1:numel(g),if m=="paired",q=abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;else,q=T.method==m&abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;end;if any(q),Z(ia,ig)=T.(v)(find(q,1));end;end;end
imagesc(ax,g-.8,A,Z);set(ax,'YDir','normal','CLim',[0 1],'XTick',g-.8,'YTick',A);xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'Vibration amplitude, A (mm)');title(ax,ttl);style(ax);
end

function draw_fail(ax,T)
cats=["EO","EO + amplitude","EO + phase","EO + amplitude + phase"];gg=sort(unique(T.g_truth_mm));M=zeros(numel(gg),numel(cats));for i=1:numel(gg),for j=1:numel(cats),q=T.g_truth_mm==gg(i)&T.failure_pattern==cats(j);if any(q),M(i,j)=T.fraction_of_failures(find(q,1));end;end;end
bar(ax,gg,M,'stacked');ylim(ax,[0 1]);xlabel(ax,'True clearance, g (mm)');ylabel(ax,'Fraction of Fixed failures');title(ax,'Fixed failure-mode composition');legend(ax,cats,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',6.5);style(ax);
end
function style(ax),set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',.75,'TickDir','in','Box','on','Layer','top','XGrid','off','YGrid','off');end
function tag(ax,s),text(ax,-.13,1.06,['(' s ')'],'Units','normalized','FontName','Times New Roman','FontSize',8.5);end
