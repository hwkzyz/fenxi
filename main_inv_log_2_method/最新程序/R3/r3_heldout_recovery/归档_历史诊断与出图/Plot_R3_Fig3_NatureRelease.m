function Files=Plot_R3_Fig3_NatureRelease()
%PLOT_R3_FIG3_NATURERELEASE R3 main figure: six-panel evidence chain.
root=fileparts(mfilename('fullpath')); S=load(fullfile(root,'output','18_paper_data_sync','r3_paper_result.mat'),'Result'); R=S.Result; P=R.P; g=sort(P.gMainMm); A=P.amplitudeMainMm;
fig=figure('Color','w','Units','centimeters','Position',[1 1 18.3 10.8]); t=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(t,1); draw_protocol(ax,P); tag(ax,'a');
ax=nexttile(t,2); draw_case(ax,root,P); tag(ax,'b');
ax=nexttile(t,3); draw_failure(ax,R.FailureByGap); tag(ax,'c');
ax=nexttile(t,4); draw_grid(ax,R.Summary,"fixed",g,A,'P_vib','Fixed-reference recovery fraction'); tag(ax,'d');
ax=nexttile(t,5); draw_grid(ax,R.Summary,"adaptive",g,A,'P_vib','Clearance-adaptive recovery fraction'); tag(ax,'e');
ax=nexttile(t,6); draw_grid(ax,R.Paired,"paired",g,A,'P_rescue','Paired rescue fraction'); tag(ax,'f');
colormap(fig,tealmap(256)); out=fullfile(root,'output','figures_final'); if ~exist(out,'dir'),mkdir(out);end; base=fullfile(out,'Fig3_R3_nature_release');
exportgraphics(fig,[base '.svg'],'ContentType','vector'); exportgraphics(fig,[base '.pdf'],'ContentType','vector'); exportgraphics(fig,[base '.emf'],'ContentType','vector'); exportgraphics(fig,[base '.png'],'Resolution',600); Files=struct('svg',[base '.svg'],'pdf',[base '.pdf'],'emf',[base '.emf'],'png',[base '.png']);
end

function draw_protocol(ax,P)
hold(ax,'on'); plot(ax,P.gAllMm,zeros(size(P.gAllMm)),'-','Color',[.75 .77 .77],'LineWidth',.7); scatter(ax,P.gCalMm,zeros(size(P.gCalMm)),25,[.32 .36 .36],'filled'); scatter(ax,P.gOperatingMm,zeros(size(P.gOperatingMm)),31,[.10 .43 .49],'filled'); scatter(ax,P.gReferenceMm,0,52,[.75 .25 .20],'o','LineWidth',1.0);
xlim(ax,[.15 1.55]);ylim(ax,[-.34 .34]);set(ax,'YTick',[],'XTick',.2:.2:1.4);xlabel(ax,'Clearance, g (mm)');title(ax,'Held-out validation design');
text(ax,.30,.18,'calibration family','FontSize',7.2,'Color',[.25 .28 .28]);text(ax,1.00,-.18,'held-out operating states','FontSize',7.2,'Color',[.05 .34 .39],'HorizontalAlignment','center');text(ax,.80,.18,'low-speed reference, g_{ref}=0.8 mm','FontSize',7.0,'Color',[.65 .18 .14],'HorizontalAlignment','center');
text(ax,.22,-.29,'calibration states  ->  F_{cal}(g,x)  ->  Fixed / Adaptive inverse','FontSize',6.5);text(ax,.22,-.34,'held-out raw curve  ->  high-speed truth record  (never enters F_{cal})','FontSize',6.5,'Color',[.05 .34 .39]);style(ax);
end

function draw_case(ax,root,P)
F=readtable(fullfile(root,'output','05_main_15db','main_detail_merged.csv'),'TextType','string');g0=.7;A0=.2;seed=P.mainHighSeeds(1);q=F.g_truth_mm==g0&F.A_true_mm==A0&F.high_seed==seed;fixed=F(q&F.method=="fixed",:);adapt=F(q&F.method=="adaptive",:);frot=fixed.f_true_hz/fixed.eo_true;x=linspace(0,.5,600);uTruth=A0*sin(2*pi*fixed.eo_true*x+fixed.phi_mapped_rad);uFixed=fixed.A_est_mm*sin(2*pi*(fixed.f_est_hz/frot)*x+fixed.phi_est_rad);uAdapt=adapt.A_est_mm*sin(2*pi*(adapt.f_est_hz/frot)*x+adapt.phi_est_rad);
plot(ax,x,uTruth,'k-','LineWidth',1.0,'DisplayName',sprintf('Truth (EO=%d, A=%.3f mm)',fixed.eo_true,A0));hold(ax,'on');plot(ax,x,uFixed,'--','Color',[.78 .29 .20],'LineWidth',1.0,'DisplayName',sprintf('Fixed (EO=%d, A=%.3f mm)',fixed.eo_est,fixed.A_est_mm));plot(ax,x,uAdapt,'-','Color',[.05 .43 .49],'LineWidth',1.0,'DisplayName',sprintf('Adaptive (EO=%d, A=%.3f mm)',adapt.eo_est,adapt.A_est_mm));xlabel(ax,'Rotor revolution');ylabel(ax,'Blade displacement, u (mm)');title(ax,'Representative paired reconstruction');ylim(ax,[-.55 .55]);legend(ax,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',5.8);style(ax);
end

function draw_grid(ax,T,m,g,A,v,ttl)
Z=NaN(numel(A),numel(g));for ia=1:numel(A),for ig=1:numel(g),if m=="paired",q=abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;else,q=T.method==m&abs(T.g_truth_mm-g(ig))<1e-12&abs(T.A_true_mm-A(ia))<1e-12;end;if any(q),Z(ia,ig)=T.(v)(find(q,1));end;end;end
imagesc(ax,1:numel(g),A,Z);set(ax,'YDir','normal','CLim',[0 1],'XTick',1:numel(g),'XTickLabel',compose('%+.1f',g-.8),'YTick',A);xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'Vibration amplitude, A (mm)');title(ax,ttl);hold(ax,'on');for xx=.5:1:(numel(g)+.5),xline(ax,xx,'Color',[.94 .95 .95],'LineWidth',.35);end;for yy=[A(1)-.025,(A(1:end-1)+A(2:end))/2,A(end)+.025],yline(ax,yy,'Color',[.94 .95 .95],'LineWidth',.35);end;style(ax);
end

function draw_failure(ax,T)
cats=["EO","EO+A","EO+phi","EO+A+phi"];gg=sort(unique(T.g_truth_mm));M=zeros(numel(gg),numel(cats));for i=1:numel(gg),for j=1:numel(cats),long=["EO","EO + amplitude","EO + phase","EO + amplitude + phase"];q=T.g_truth_mm==gg(i)&T.failure_pattern==long(j);if any(q),M(i,j)=T.fraction_of_failures(find(q,1));end;end;end
h=bar(ax,gg-.8,M,'stacked','BarWidth',.72);C=[.29 .40 .47;.78 .36 .20;.79 .63 .25;.46 .35 .50];for j=1:numel(h),h(j).FaceColor=C(j,:);h(j).EdgeColor='none';end;ylim(ax,[0 1]);xlabel(ax,'Clearance mismatch, \Delta g (mm)');ylabel(ax,'Fraction of Fixed failures');title(ax,'Fixed failure composition');legend(ax,h,cats,'Location','northeast','Box','off','FontSize',5.8);style(ax);
end
function style(ax),set(ax,'FontName','Arial','FontSize',7.2,'LineWidth',.6,'TickDir','in','Box','on','Layer','top','XGrid','off','YGrid','off');end
function tag(ax,s),text(ax,-.12,1.05,['(' s ')'],'Units','normalized','FontName','Arial','FontSize',8.5,'FontWeight','bold');end
function c=tealmap(n),x=[0 .25 .5 .75 1];a=[.95 .96 .95;.82 .90 .89;.60 .78 .78;.25 .57 .59;.03 .32 .36];c=interp1(x,a,linspace(0,1,n),'pchip');end
