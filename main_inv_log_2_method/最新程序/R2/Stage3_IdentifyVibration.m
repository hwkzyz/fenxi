function estimate = Stage3_IdentifyVibration(caseID,forceRegenerate)
%STAGE3_IDENTIFYVIBRATION 只执行 R2 振动/EO/间隙盲辨识。
if nargin<1, caseID=1; end; if nargin<2, forceRegenerate=false; end
here=fileparts(mfilename('fullpath')); root=fileparts(here); addpath(fullfile(root,'local_func'),'-begin');
outDir=fullfile(here,'r2_profiled_recovery','output','stages'); if ~exist(outDir,'dir'), mkdir(outDir); end
file=fullfile(outDir,sprintf('R2_identification_case_%02d.mat',caseID));
if exist(file,'file') && ~forceRegenerate, S=load(file,'estimate'); estimate=S.estimate; return; end
% 辨识阶段只读取已经生成的波形和响应库，不再隐式调用前两个阶段。
waveFile=fullfile(outDir,sprintf('R2_waveform_case_%02d.mat',caseID));
libFile=fullfile(outDir,'R2_gap_library.mat');
if ~exist(waveFile,'file'), error('缺少波形文件，请先运行 Stage2_GenerateWaveform(%d)。',caseID); end
if ~exist(libFile,'file'), error('缺少间隙库，请先运行 Stage1_BuildGapLibrary。'); end
sw=load(waveFile,'waveform'); waveform=sw.waveform; m=waveform.mapped;
ctx=load_inv_log_2_project_context(root);
bundle=struct('t',m.t_v,'x',m.x_v,'V',m.V_a,'theta',m.theta_v,'sensorId',m.S_v,'turnId',m.rev_v,'Wvp',ones(size(m.V_a)),'Wfull',ones(size(m.V_a)));
response=struct('evalF',@(g,x)eval_gap_template(waveform.calibration.templateModel,g,x),'evalFx',@(g,x)eval_gap_derivative(waveform.calibration.templateModel,g,x));
gb=[min(ctx.gapList),max(ctx.gapList)]; bc=struct('gapBounds',gb,'gapGrid',gb(1):.02:gb(2),'eoGrid',1:30,'topK',5,'rotHz',waveform.cfg.RPM_high/60,'minGradient',1e-6,'minSamples',30,'dxBound',.2,'ampCoeffBound',1.5,'deltaFBound',2,'maxIter',250,'ambiguityMargin',.01);
B=run_unified_blind_dynamic_backend(bundle,response,bc); estimate=struct('caseID',caseID,'backend',B,'truth',waveform.truth); save(file,'estimate','-v7.3');
end
