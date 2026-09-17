function Files=Plot_R3_Fig3_NatureFinal()
%PLOT_R3_FIG3_NATUREFINAL Final 2x3 evidence-chain figure for R3.
root=fileparts(mfilename('fullpath')); S=load(fullfile(root,'output','18_paper_data_sync','r3_paper_result.mat'),'Result'); R=S.Result; P=R.P; g=sort(P.gMainMm); A=P.amplitudeMainMm;
fig=figure('Color','w','Units','centimeters','Position',[1 1 17 10.8]); t=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(t,1); protocol(ax,P); tag(ax,'a');
ax=nexttile(t,2); representative(ax,root,P); tag(ax,'b');
ax=nexttile(t,3); heat(ax,R.Summary,"fixed",g,A,'P_vib','Fixed reference'); tag(ax,'c');
ax=nexttile(t,4); heat(ax,R.Summary,"adaptive",g,A,'P_vib','Clearance adaptive'); tag(ax,'d');
ax=nexttile(t,5); heat(ax,R.Paired,"paired",g,A,'P_rescue','Paired adaptive rescue'); tag(ax,'e');
ax=nexttile(t,6); failure(ax,R.FailureByGap); tag(ax,'f');
colormap(fig,tealmap(256)); title(t,'Synchronous-vibration recovery under clearance mismatch','FontName','Arial','FontSize',9,'FontWeight','normal');
out=fullfile(root,'output','figures_final'); if ~exist(out,'dir'),mkdir(out);end; base=fullfile(out,'Fig3_R3_nature_final');
exportgraphics(fig,[base '.png'],'Resolution',600); exportgraphics(fig,[base '.pdf'],'ContentType','vector'); exportgraphics(fig,[base '.emf'],'ContentType','vector'); Files=struct('png',[base '.png'],'pdf',[base '.pdf'],'emf',[base '.emf']);
end

function protocol(ax,P)
hold(ax,'on'); plot(ax,P.gAllMm,zeros(size(P.gAllMm)),'-','Color',[.75 .77 .77],'LineWidth',.8); scatter(ax,P.gCalMm,zeros(size(P.gCalMm)),28,[.32 .36 .36],'filled'); scatter(ax,P.gOperatingMm,zeros(size(P.gOperatingMm)),34,[.10 .43 .49],'filled'); scatter(ax,P.gReferenceMm,0,55,[.75 .25 .20],'o','LineWidth',1.1);
xlim(ax,[.15 1.55]);ylim(ax,[-.28 .28]);set(ax,'YTick',[],'XTick',.2:.2:1.4);xlabel(ax,'Clearance, g (mm)');title(ax,'Frozen calibration and held-out test design');text(ax,.31,.17,'calibration family','FontSize',7.5,'Color',[.25 .28 .28]);text(ax,1.00,-.16,'held-out operating states','FontSize',7.5,'Color',[.05 .34 .39],'HorizontalAlignment','center');text(ax,.8,.18,'g_{ref}=0.8 mm','FontSize',7.5,'Color',[.65 .18 .14],'HorizontalAlignment','center');text(ax,.22,-.25,'raw held-out curve  \rightarrow  noisy record  \rightarrow  calibration-only inverse','FontSize',6.8);style(ax);
end

function representative(ax,root,P)
F=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');g0=.7;A0=.2;seed=P.mainHighSeeds(1);q=F.g_truth_mm==g0&F.A_true_mm==A0&F.high_seed==seed;fixed=F(q&F.method=="fixed",:);adapt=F(q&F.method=="adaptive",:);
x=linspace(0,1,600);phi=P.phaseValuesRad(1);uTruth=A0*sin(2*pi*5*x+phi);uFixed=fixed.A_est_mm*sin(2*pi*(fixed.eo_est/10)*5*x+fixed.phi_est_rad);uAdapt=adapt.A_est_mm*sin(2*pi*(adapt.eo_est/10)*5*x+adapt.phi_est_rad);
plot(ax,x,uTruth,'k-','LineWidth',1.0,'DisplayName','Truth');hold(ax,'on');plot(ax,x,uFixed,'--','Color',[.78 .29 .20],'LineWidth',1.0,'DisplayName','Fixed reconstruction');plot(ax,x,uAdapt,'-','Color',[.05 .43 .49],'LineWidth',1.0,'DisplayName','Adaptive reconstruction');xlabel(ax,'Rotor revolution, \theta/(2\pi)');ylabel(ax,'Synchronous displacement, u (mm)');title(ax,'Representative paired case: g=0.7 mm, A=0.20 mm');ylim(ax,[-.6 .6]);legend(ax,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',6.5);text(ax,.03,.94,sprintf('Truth: EO10, A=%.2f mm\\nFixed: EO%d, A=%.3f mm\\nAdaptive: EO%d, A=%.3f mm',A0,fixed.eo_est,fixed.A_est_mm,adapt.eo_est,adapt.A_est_mm),'Units','normalized','VerticalAlignment','top','FontSize',6.7,'BackgroundColor','w','Margin',2);style(ax);
end

function heat(ax,T,m,g,A,v,ttl)
Z=NaN(numel(A),numel(g));for ia=1:numel(A),for ig=1:numel(g),if m=="paired",q=abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;else,q=T.method==m&abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;end;if any(q),Z(ia,ig)=T.(v)(find(q,1));end;end;end
imagesc(ax,g-.8,A,Z);set(ax,'YDir','normal','CLim',[0 1],'XTick',g-.8,'YTick',A);xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'Vibration amplitude, A (mm)');title(ax,ttl);hold(ax,'on');for x=g-.8,xline(ax,x,'Color','w','LineWidth',.6);end;for y=A,yline(ax,y,'Color','w','LineWidth',.6);end;style(ax);
end

function failure(ax,T)
cats=["EO","EO + amplitude","EO + phase","EO + amplitude + phase"];gg=sort(unique(T.g_truth_mm));M=zeros(numel(gg),numel(cats));for i=1:numel(gg),for j=1:numel(cats),q=T.g_truth_mm==gg(i)&T.failure_pattern==cats(j);if any(q),M(i,j)=T.fraction_of_failures(find(q,1));end;end;end
bar(ax,gg-.8,M,'stacked','BarWidth',.75);ylim(ax,[0 1]);xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'Fraction of Fixed failures');title(ax,'Fixed failure modes');legend(ax,cats,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',5.8);colors=[.27 .37 .44;.78 .36 .20;.79 .63 .25;.42 .30 .47];h=findobj(ax,'Type','Bar');for i=1:min(numel(h),size(colors,1)),h(i).FaceColor=colors(size(colors,1)-i+1,:);end;style(ax);
end
function style(ax),set(ax,'FontName','Arial','FontSize',7.5,'LineWidth',.65,'TickDir','in','Box','on','Layer','top','XGrid','off','YGrid','off');end
function tag(ax,s),text(ax,-.12,1.04,['(' s ')'],'Units','normalized','FontName','Arial','FontSize',8.5,'FontWeight','bold');end
function c=tealmap(n),x=[0 .5 1];a=[.93 .95 .94;.43 .69 .70;.03 .32 .36];c=interp1(x,a,linspace(0,1,n),'pchip');end
