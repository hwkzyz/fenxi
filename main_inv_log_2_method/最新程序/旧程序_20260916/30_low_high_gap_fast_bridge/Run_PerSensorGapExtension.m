function Result = Run_PerSensorGapExtension()
%RUN_PERSENSORGAPEXTENSION Smoke test for the optional per-sensor extension.

root=fileparts(mfilename('fullpath')); mainDir=fileparts(root);
addpath(mainDir,'-begin'); addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);
trust=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN, ...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),trust);
base=make_inv_log_2_demo_case(ctx,.2); c=base.cfgCase;
c.RPM_low=min(c.RPM_high,300); c.NumRevs_low=8; c.snrDb=20;
c.route30BaseSamples=600; c.route30ForwardModel="absolute";
c.perSensorFrequencySeed=[700,1200];
gLow=.8; gShared=.5; gSensor=[.40,.50,.60];
lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z), ...
    c,lib.domain,20,20260807);
low=aggregate_low_speed_template(lowData,c.alpha_k,c.R_tip,lib.domain,lib.xGrid);
lowState=estimate_highspeed_static_gap_raw(low.mapped,lib,c);

c.A_true=[.25,.15]; c.f_true=[700,1200]; c.phi_true=[pi/4,-pi/3];
d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gShared,z), ...
    c.RPM_high,c.NumRevs_high,c.fs,c.R_tip,c.alpha_k,0,lib.domain, ...
    c.A_true,c.f_true,c.phi_true,'noise_ratio');
sig=rms(d.V_clean-min(d.V_clean)); rng(20260807,'twister');
d.V_cap=d.V_clean+sig/10^(c.snrDb/20)*randn(size(d.V_clean));
m=map_highspeed_to_space(d,c.alpha_k,c.R_tip,lib.domain,.02);
u=c.A_true(1)*sin(2*pi*c.f_true(1)*m.t_v+c.phi_true(1))+ ...
    c.A_true(2)*sin(2*pi*c.f_true(2)*m.t_v+c.phi_true(2));
baseClean=eval_gap_template(lib,gShared,m.x_v-u);
noise=m.V_a-baseClean;
sensorValues=unique(m.S_v,'stable');
for s=1:numel(sensorValues)
    idx=m.S_v==sensorValues(s);
    m.V_a(idx)=eval_gap_template(lib,gSensor(s),m.x_v(idx)-u(idx))+noise(idx);
end

fit=run_per_sensor_gap_extension(m,lib,c,lowState);
Result=struct('fit',fit,'g_true',gSensor,'outputDir',fullfile(root,'output','per_sensor_extension'));
out=Result.outputDir; if ~exist(out,'dir'),mkdir(out);end
T=table(sensorValues(:),gSensor(:),fit.sensor_gap_used(:), ...
    fit.sensor_gap_increment(:),'VariableNames',{'sensor','g_true','g_est','dg_est'});
writetable(T,fullfile(out,'per_sensor_gap_extension.csv')); save(fullfile(out,'per_sensor_gap_extension.mat'),'T','fit');
disp(T);
fprintf('f_est=%g,%g rmse=%g\n',fit.f_id(1),fit.f_id(2),fit.rmse);
if max(abs(T.g_est-T.g_true))>.03 || max(abs(fit.f_id-c.f_true))>2
    error('Per-sensor extension smoke test failed.');
end
end
