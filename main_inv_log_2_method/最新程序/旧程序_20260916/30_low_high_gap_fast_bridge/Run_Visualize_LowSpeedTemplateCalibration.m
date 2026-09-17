function Out = Run_Visualize_LowSpeedTemplateCalibration(snrDb,numRevs,sensorId)
%RUN_VISUALIZE_LOWSPEEDTEMPLATECALIBRATION Visualize low-speed template marking.
if nargin<1 || isempty(snrDb),snrDb=15;end
if nargin<2 || isempty(numRevs),numRevs=20;end
if nargin<3 || isempty(sensorId),sensorId=1;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
cfg=make_inv_log_2_demo_case(ctx,.2);cfg=cfg.cfgCase;cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=numRevs;
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,.5,z),cfg,lib.domain,snrDb,941001);
low=build_low_speed_templates_step04(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',9,'minCoverage',.90));
mapped=low.mapped;idxS=mapped.S_v==sensorId;x=mapped.x_v(idxS);v=mapped.V_a(idxS);xg=lib.xGrid(:);tpl=low.templateBySensor(:,ismember(low.sensorIds,sensorId));truth=eval_gap_template(lib,.5,xg);d=median(diff(xg));
rev=unique(mapped.rev_v(idxS),'stable');Y=nan(numel(xg),numel(rev));
for ir=1:numel(rev)
    ii=idxS & mapped.rev_v==rev(ir);xx=mapped.x_v(ii);vv=mapped.V_a(ii);[xx,ord]=sort(xx);vv=vv(ord);[xx,keep]=unique(xx,'stable');vv=vv(keep);
    if numel(xx)>3,Y(:,ir)=interp1(xx,vv,xg,'pchip',NaN);end
end
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 12.8]);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on;box on;tt=lowData.t;plot(tt,lowData.V_cap,'Color',[.70 .70 .70],'LineWidth',.45);xlim([tt(1),tt(min(end,round(numel(tt)*.18)))]);xlabel('Time (s)');ylabel('Voltage (V)');title('(a) Noisy low-speed record');
nexttile;hold on;box on;scatter(x,v,3,mapped.rev_v(idxS),'.');xlabel('Spatial coordinate x (mm)');ylabel('Voltage (V)');title('(b) OPR-aligned point cloud');colormap(turbo);cb=colorbar;cb.Label.String='Turn index';
nexttile;hold on;box on;plot(xg,Y,'Color',[.78 .78 .78],'LineWidth',.35);plot(xg,tpl,'k','LineWidth',1.4);plot(xg,truth,'r--','LineWidth',1.0);xlabel('Spatial coordinate x (mm)');ylabel('Voltage (V)');title('(c) Turn aggregation and template');legend({'Per-turn interpolation','Calibrated template','Noise-free truth'},'Location','best');
nexttile;hold on;box on;plot(xg,gradient(tpl,d),'k','LineWidth',1.3);plot(xg,gradient(truth,d),'r--','LineWidth',1.0);xlabel('Spatial coordinate x (mm)');ylabel('dV/dx (V/mm)');title('(d) Template derivative');legend({'Estimated','Truth'},'Location','best');
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');set(findall(fig,'-property','FontSize'),'FontSize',8);set(findall(fig,'Type','axes'),'TickDir','in','FontSize',8);
outDir=fullfile(root,'output','low_speed_template_visualization');if ~exist(outDir,'dir'),mkdir(outDir);end
tag=sprintf('snr_%g_revs_%d_sensor_%d',snrDb,numRevs,sensorId);png=fullfile(outDir,[tag,'.png']);exportgraphics(fig,png,'Resolution',240);try,exportgraphics(fig,fullfile(outDir,[tag,'.pdf']),'ContentType','vector');end
Out=struct('snrDb',snrDb,'numRevs',numRevs,'sensorId',sensorId,'figure',png,'templateRmse_V',sqrt(mean((tpl-truth).^2)),'derivativeRmse_V_per_mm',sqrt(mean((gradient(tpl,d)-gradient(truth,d)).^2)));disp(Out);
end
