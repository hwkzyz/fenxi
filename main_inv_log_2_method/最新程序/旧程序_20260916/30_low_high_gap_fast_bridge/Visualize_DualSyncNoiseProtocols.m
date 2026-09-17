function T=Visualize_DualSyncNoiseProtocols()
%VISUALIZE_DUALSYNCNOISEPROTOCOLS Compare 15 dB noise definitions visually.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','noise_protocol_waveforms');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_high=8;cfg.NumRevs_low=20;
cfg.route30ForwardModel="low_increment";cfg.A_true=[.15 .10];
fRot=cfg.RPM_high/60;cfg.f_true=[10 26]*fRot;cfg.phi_true=[pi/4 -pi/3];
gLow=.5;gHigh=.6;path=make_truth_model(lib,cfg,gLow);

% Weak two-synchronous condition and its no-vibration counterpart.
high=simulate_highspeed_from_low_increment(path,gHigh-gLow,cfg,Inf,'fixed_std',1);
highMap=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
highStatic=high;highStatic.V_cap=high.V_static;
staticMap=map_highspeed_to_space(highStatic,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
vClean=highMap.V_a(:);vStatic=staticMap.V_a(:);vVib=vClean-vStatic;

% All protocols use the same unit-noise realization, isolating only sigma.
rng(424242,'twister');z=randn(size(vClean));snrDb=15;
sFull=rms(vClean-mean(vClean));sVib=rms(vVib);
cfgRef=cfg;cfgRef.A_true=[.25 .15];cfgRef.f_true=[10 26]*fRot;cfgRef.phi_true=[pi/4 -pi/3];
refHigh=simulate_highspeed_from_low_increment(path,.1,cfgRef,Inf,'fixed_std',1);
refMap=map_highspeed_to_space(refHigh,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
refStatic=refHigh;refStatic.V_cap=refHigh.V_static;
refStaticMap=map_highspeed_to_space(refStatic,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
sReferenceVib=rms(refMap.V_a-refStaticMap.V_a);
sStatic=rms(eval_gap_template(lib,.5,lib.xGrid)-mean(eval_gap_template(lib,.5,lib.xGrid)));

labels=["Fixed reference vibration";"Full waveform SNR";"Vibration-increment SNR";"Static waveform SNR"];
sigmas=[sReferenceVib;sFull;sVib;sStatic]/10^(snrDb/20);
actualFull=20*log10(sFull./sigmas);actualVib=20*log10(sVib./sigmas);
T=table(labels,sigmas,actualFull,actualVib,'VariableNames',...
    {'noise_definition','noise_std_V','full_waveform_snr_db','vibration_increment_snr_db'});
writetable(T,fullfile(outDir,'weak_dualsync_noise_protocol_levels_15dB.csv'));

% Show a single sensor passage.  It is sufficient to display local waveform
% deformation and avoids hiding noise in the long multi-revolution record.
sel=highMap.S_v==1 & highMap.rev_v==1;x=highMap.x_v(sel);[x,ord]=sort(x);
v0=vClean(sel);v0=v0(ord);vs=vStatic(sel);vs=vs(ord);dv=v0-vs;zz=z(sel);zz=zz(ord);
f=figure('Color','w','Units','centimeters','Position',[1 1 17 12.5]);
tl=tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for i=1:4
    vn=v0+sigmas(i)*zz;dn=dv+sigmas(i)*zz;
    nexttile(i);plot(x,v0,'k-','LineWidth',1.0);hold on;plot(x,vn,'Color',[.15 .45 .75],'LineWidth',.75);
    title(sprintf('(%c) %s',char('a'+i-1),labels(i)),'FontWeight','normal');
    if i==1,ylabel('Voltage (V)');end
    if i==1,legend({'Clean','Noisy'},'Location','best');end
    nexttile(i+4);plot(x,dv,'k-','LineWidth',1.0);hold on;plot(x,dn,'Color',[.75 .25 .20],'LineWidth',.75);
    if i==1,ylabel('Vibration increment (V)');end
    xlabel('Spatial coordinate (mm)');
end
ax=findall(f,'Type','axes');for i=1:numel(ax)
    set(ax(i),'FontName','Times New Roman','FontSize',8,'TickDir','in','Box','on');
end
sgtitle(tl,sprintf('Weak dual-sync: A = [0.15, 0.10] mm; nominal 15 dB; same noise realization'),...
    'FontName','Times New Roman','FontSize',9,'FontWeight','normal');
exportgraphics(f,fullfile(outDir,'weak_dualsync_noise_protocols_15dB.png'),'Resolution',300);
exportgraphics(f,fullfile(outDir,'weak_dualsync_noise_protocols_15dB.pdf'),'ContentType','vector');
try,exportgraphics(f,fullfile(outDir,'weak_dualsync_noise_protocols_15dB.emf'),'ContentType','vector');catch,end
close(f);disp(T);
end

function model=make_truth_model(lib,cfg,gLow)
T=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg.alpha_k));
cal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,...
    'sensorIds',1:numel(cfg.alpha_k));model=build_path_template_model(lib,T,cal);
end
