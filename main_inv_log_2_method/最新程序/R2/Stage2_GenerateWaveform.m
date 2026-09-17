function waveform = Stage2_GenerateWaveform(caseID,forceRegenerate)
%STAGE2_GENERATEWAVEFORM 生成并缓存一个 R2 高低速波形案例。
if nargin < 1, caseID = 1; end
if nargin < 2, forceRegenerate = false; end
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'local_func'),'-begin'); addpath(root,'-begin');
outDir = fullfile(here,'r2_profiled_recovery','output','stages');
if ~exist(outDir,'dir'), mkdir(outDir); end
file = fullfile(outDir,sprintf('R2_waveform_case_%02d.mat',caseID));
if exist(file,'file') && ~forceRegenerate, S=load(file,'waveform'); waveform=S.waveform; return; end
libFile=fullfile(outDir,'R2_gap_library.mat');
if ~exist(libFile,'file'), error('缺少间隙库，请先运行 Stage1_BuildGapLibrary。'); end
sl=load(libFile,'lib'); lib=sl.lib;
ctx = load_inv_log_2_project_context(root);
% 这里只需要分析配置；不要调用 make_inv_log_2_demo_case，那个函数还会
% 额外生成一套演示高速波形，属于旧的“一步式”流程。
cfg = ctx.cfgAna;
cfg.g_holdout = 0.2;
cfg.A_true = 0.2;
cfg.f_true = 500;
cfg.phi_true = pi/4;
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_low=20; cfg.NumRevs_high=8;
cfg.route30ForwardModel="low_increment"; cfg.route30FrequencyStructureMode="single_sync";
eoList=[7 10 13 16 19 26]; ampList=[.20 .25 .16 .22 .18 .24];
if caseID>numel(eoList), error('caseID 超出 1:%d',numel(eoList)); end
gLow=.5; gHigh=.6; truthT=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg.alpha_k));
truthCal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));
truthModel=build_path_template_model(lib,truthT,truthCal);
low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),cfg,lib.domain,Inf,88000+caseID,'fixed_std');
calibration=calibrate_inv_log_2_low_speed(low,lib,cfg);
cfg.A_true=ampList(caseID); cfg.f_true=eoList(caseID)*cfg.RPM_high/60; cfg.phi_true=pi/4;
high=simulate_highspeed_from_low_increment(truthModel,gHigh-gLow,cfg,Inf,'fixed_std',89000+caseID);
mapped=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02,.30);
waveform=struct('caseID',caseID,'low',low,'high',high,'mapped',mapped,'calibration',calibration, ...
 'cfg',cfg,'truth',struct('eo',eoList(caseID),'A',ampList(caseID),'gLow',gLow,'gHigh',gHigh));
save(file,'waveform','-v7.3'); fprintf('R2 波形已保存：%s\n',file);
end
