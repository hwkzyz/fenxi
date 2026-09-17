function Result=Run_StructuredGapModelAblation(snrDb,seedOffset,truthMode,caseFilter,methodFilter)
%RUN_STRUCTUREDGAPMODELABLATION Compare gap-model choices without changing solvers.
% Methods: fixed_low_gap, absolute_joint, low_template_increment.
if nargin<1,snrDb=15;end
if nargin<2,seedOffset=0;end
if nargin<3,truthMode="low_increment";end
if nargin<4,caseFilter="all";end
if nargin<5,methodFilter="all";end
truthMode=string(truthMode);caseFilter=string(caseFilter);methodFilter=string(methodFilter);
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=configure_structured_identification(base.cfgCase,"experiment_current");
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.NumRevs_high=8;cfg.snrDb=snrDb;
fRot=cfg.RPM_high/60;cases=ablation_cases(fRot);
if ~any(caseFilter=="all")
    cases=cases(ismember([cases.id],caseFilter));
end
if isempty(cases),error('No ablation case matches caseFilter.');end
allMethodName=["fixed_low_gap","absolute_joint","low_template_increment"];
methodMask=methodFilter=="all"|ismember(allMethodName,methodFilter);
if ~any(methodMask),error('No ablation method matches methodFilter.');end
rows=repmat(row0(),numel(cases)*nnz(methodMask),1);at=0;
runTag=sprintf('structured_gap_ablation_%gdb_seed%d_%s',snrDb,seedOffset,truthMode);
if ~any(caseFilter=="all"),runTag=runTag+"_case_"+strjoin(caseFilter,"-");end
if ~any(methodFilter=="all"),runTag=runTag+"_method_"+strjoin(methodFilter,"-");end
out=fullfile(root,'output',runTag);if ~exist(out,'dir'),mkdir(out);end
for ic=1:numel(cases)
    c=cases(ic);lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,.2,z),...
        cfg,lib.domain,snrDb,240000+seedOffset+ic);
    low=aggregate_low_speed_template(lowData,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid);
    needsPathModel=truthMode=="low_increment"||methodMask(3);
    if needsPathModel
        pathCal=calibrate_low_speed_path(lib,low,struct('xDomain',lib.domain,...
            'fitMask',[true false false false],'fitGainOffset',false));
        pathModel=build_path_template_model(lib,low.templateLow,pathCal);
    else
        pathCal=struct('g0',.2);pathModel=[];
    end
    gHigh=c.gHigh;delta=gHigh-pathCal.g0;
    % The path-increment simulator reads the truth waveform from cfg. Keep
    % these fields synchronized with the case rather than inheriting the
    % demo case's unrelated defaults.
    cfgCase=cfg;cfgCase.A_true=c.A;cfgCase.f_true=c.f;cfgCase.phi_true=c.phi;
    if truthMode=="low_increment"
        highInc=simulate_highspeed_from_low_increment(pathModel,delta,cfgCase,snrDb,'snr_db',250000+seedOffset+ic);
    else
        highInc=[];
    end
    if truthMode=="absolute"
        high=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
            cfg.RPM_high,8,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,c.A,c.f,c.phi,'noise_ratio');
        sig=rms(high.V_clean-min(high.V_clean));rng(250100+seedOffset+ic,'twister');
        high.V_cap=high.V_clean+sig/10^(snrDb/20)*randn(size(high.V_clean));
    else
        high=highInc;
    end
    map=map_highspeed_to_space(high,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    methodCfg={fixed_gap_cfg(cfg),cfg,cfg};
    methodCfg{1}.route30ForwardModel="absolute";methodCfg{2}.route30ForwardModel="absolute";
    methodCfg{3}.route30ForwardModel="low_increment";methodLib={lib,lib,pathModel};
    methodName=allMethodName;
    for im=1:3
        if ~methodMask(im),continue;end
        at=at+1;ticFit=tic;
        fit=run_inv_log_2_structured_main(low,map,methodLib{im},methodCfg{im},c.mode);
        elapsed=toc(ticFit);
        [ferr,ampErr,phaseErr]=parameter_errors(fit,c,fRot);
        rows(at)=struct('case_id',c.id,'mode',c.mode,'method',methodName(im),...
            'snr_db',snrDb,'truth_mode',truthMode,'frequency_error_hz',ferr,...
            'gap_error_mm',abs(fit.g_used-gHigh),'amplitude_true',string(mat2str(c.A,7)),...
            'amplitude_est',string(mat2str(fit.A_id,7)),...
            'max_amplitude_relative_error',ampErr,'max_phase_error_rad',phaseErr,'rmse',fit.rmse,...
            'identification_status',fit.identification_status,'elapsed_s',elapsed,...
            'success',ferr<=1&&abs(fit.g_used-gHigh)<=.05&&ampErr<=.10);
        % Long mixed-mode fits are checkpointed so an interrupted matrix
        % run does not discard already completed ablation rows.
        Audit=struct2table(rows(1:at));
        writetable(Audit,fullfile(out,'structured_gap_ablation.csv'));
        save(fullfile(out,'structured_gap_ablation.mat'),'Audit');
    end
end
Audit=struct2table(rows(1:at));disp(Audit);
Result=struct('Audit',Audit,'outputDir',out);
end

function c=fixed_gap_cfg(c)
c.route30GapHalfWidthMm=0;c.route30GapCount=1;c.mainVpGapHalfWidth=0;c.mainVpGapN=1;
c.route30SingleGapHalfWidthMm=0;c.dualSyncGapCount=1;c.syncAsyncGapCount=1;
c.asyncSingleDirectGapCount=1;
end

function [ferr,ampErr,phaseErr]=parameter_errors(fit,c,fRot)
[fTrue,it]=sort(c.f(:).');[fEst,ie]=sort(fit.f_id(:).');
if numel(fTrue)~=numel(fEst)
    ferr=Inf;ampErr=Inf;phaseErr=Inf;return;
end
Atrue=c.A(:).';Atrue=Atrue(it);Aest=fit.A_id(:).';Aest=Aest(ie);
pTrue=c.phi(:).';
% Synchronous solvers use shaft angle referenced to the OPR pulse, while
% the simulator's truth phase is referenced to t=0. The simulated OPR
% delay is 0.08 revolution, so convert synchronous truth phases before
% reporting phase error. Asynchronous components remain on the time basis.
syncMask=false(size(pTrue));
if c.mode=="single_sync"||c.mode=="dual_sync_sync",syncMask(:)=true;
elseif c.mode=="dual_sync_async",syncMask(1)=true;
end
pTrue(syncMask)=pTrue(syncMask)+2*pi*(c.f(syncMask)/fRot)*0.08;
pTrue=pTrue(it);pEst=fit.phi_id(:).';pEst=pEst(ie);
ferr=max(abs(fEst-fTrue));
ampErr=max(abs(Aest-Atrue)./max(abs(Atrue),eps));
phaseErr=max(abs(atan2(sin(pEst-pTrue),cos(pEst-pTrue))));
end

function C=ablation_cases(fRot)
p=[pi/4 -pi/3];C=[...
 struct('id',"single_sync_500",'mode',"single_sync",'gHigh',.8,'A',.1,'f',500,'phi',p(1)),...
 struct('id',"single_async_733",'mode',"single_async",'gHigh',.8,'A',.15,'f',733,'phi',p(1)),...
 struct('id',"dual_sync_500_1300",'mode',"dual_sync_sync",'gHigh',.8,'A',[.25 .15],'f',[500 1300],'phi',p),...
 struct('id',"sync_async_500_733",'mode',"dual_sync_async",'gHigh',.8,'A',[.25 .15],'f',[10*fRot 733],'phi',p)];
end

function r=row0()
r=struct('case_id',"",'mode',"",'method',"",'snr_db',NaN,'truth_mode',"",...
    'frequency_error_hz',NaN,'gap_error_mm',NaN,'amplitude_true',"",'amplitude_est',"",...
    'max_amplitude_relative_error',NaN,'max_phase_error_rad',NaN,'rmse',NaN,...
    'identification_status',"",'elapsed_s',NaN,'success',false);
end
