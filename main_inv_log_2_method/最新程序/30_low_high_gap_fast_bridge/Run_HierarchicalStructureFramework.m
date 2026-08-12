function [Detail,Summary,Confusion,ReferenceNoise]=Run_HierarchicalStructureFramework(numSeeds,snrList,layers,caseModes)
%RUN_HIERARCHICALSTRUCTUREFRAMEWORK Two-layer paper simulation entry.
% known: structure-informed sub-solver validation; blind: unknown-structure
% hierarchical routing.  Both layers share exactly the same low/high data.
if nargin<1||isempty(numSeeds),numSeeds=1;end
if nargin<2||isempty(snrList),snrList=25;end
if nargin<3||isempty(layers),layers=["known","blind"];end
if nargin<4||isempty(caseModes)
    caseModes=["single_sync","single_async","dual_sync_sync","dual_sync_async"];
end
layers=string(layers);caseModes=string(caseModes);
frameworkVersion="v3_declared_frequency_domain";

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','hierarchical_structure_framework');
if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="low_increment";cfg0.route30LowTemplateMethod="adaptive_sg";
cfg0.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm',...
    [.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
% This is a declared analysis domain, shared by every truth case and never
% changed from the generated frequencies.  The solver still does not know
% the actual EO/frequency inside this interval.
cfg0.route30SingleFrequencyRangeHz=[500 1300];cfg0.route30DualFrequencyRangeHz=[500 1300];
cfg0.route30GapHalfWidthMm=.20;cfg0.structuredComponentAmplitudeFloorMm=.05;
cfg0.structuredRequireMinRevolutions=true;cfg0.structuredMinRevolutions=8;
cfg0.dualSyncSearchStrategy="coarse_top_m_multistart";
cfg0.dualSyncGapCount=11;cfg0.dualSyncKeepPerGap=80;
cfg0.dualSyncIteratedReplayCount=2100;cfg0.dualSyncRefineCount=50;
cfg0.dualSyncTopPairCount=40;cfg0.dualSyncStartsPerTopPair=2;

gLow=.50;deltaTrue=.10;gHigh=gLow+deltaTrue;fRot=cfg0.RPM_high/60;
truthT=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg0.alpha_k));
truthCal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,...
    'sensorIds',1:numel(cfg0.alpha_k));
truthModel=build_path_template_model(lib,truthT,truthCal);
referenceCfg=cfg0;referenceCfg.A_true=[.25 .15];referenceCfg.f_true=[10 26]*fRot;
referenceCfg.phi_true=[pi/4 -pi/3];
sigmaRef=reference_noise_std(truthModel,deltaTrue,referenceCfg,lib.domain,snrList);
ReferenceNoise=table(snrList(:),sigmaRef(:),...
    'VariableNames',{'reference_snr_db','fixed_noise_std_V'});
writetable(ReferenceNoise,fullfile(outDir,'reference_noise.csv'));

nMax=numSeeds*numel(snrList)*numel(caseModes)*numel(layers);
rows=repmat(empty_row(),nMax,1);ir=0;checkpoint=fullfile(outDir,'detail_checkpoint.csv');
if exist(checkpoint,'file')
    old=readtable(checkpoint,'TextType','string');
    if ismember('framework_version',old.Properties.VariableNames)
        keep=old.framework_version==frameworkVersion&ismember(old.seed,1:numSeeds)&...
            ismember(old.reference_snr_db,snrList)&ismember(old.truth_mode,caseModes)&...
            ismember(old.validation_layer,layers);
        oldRows=table2struct(old(keep,:));nOld=min(numel(oldRows),nMax);
        rows(1:nOld)=oldRows(1:nOld);ir=nOld;fprintf('Resuming from %d runs.\n',ir);
    end
end

for isnr=1:numel(snrList)
    snrDb=snrList(isnr);sigma=sigmaRef(isnr);
    for seed=1:numSeeds
        lowSeed=981000+seed;
        lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
            cfg0,lib.domain,sigma,lowSeed,'fixed_std');
        calTic=tic;C=calibrate_inv_log_2_low_speed(lowData,lib,cfg0);calTime=toc(calTic);
        for im=1:numel(caseModes)
            truthMode=caseModes(im);truth=case_truth(truthMode,fRot);
            cfgData=cfg0;cfgData.A_true=truth.A;cfgData.f_true=truth.f;cfgData.phi_true=truth.phi;
            highSeed=982000+100*im+seed;
            high=simulate_highspeed_from_low_increment(truthModel,deltaTrue,cfgData,...
                sigma,'fixed_std',highSeed);
            highMap=map_highspeed_to_space(high,cfgData.alpha_k,cfgData.R_tip,lib.domain,.02);
            for il=1:numel(layers)
                layer=layers(il);
                if is_completed(rows,ir,seed,snrDb,truthMode,layer),continue;end
                cfg=cfgData;
                if layer=="known"
                    cfg.route30FrequencyStructureMode=truthMode;
                    cfg.route30StructuredFastBudgets=~startsWith(truthMode,"dual_");
                elseif layer=="blind"
                    cfg.route30FrequencyStructureMode="unknown";
                    cfg.unknownStructureNoiseStdV=sigma;
                else
                    error('Unknown validation layer: %s.',layer);
                end
                fitTic=tic;
                try
                    fit=run_inv_log_2_high_with_calibration(C,highMap,cfg);fitTime=toc(fitTic);
                    row=pack_row(frameworkVersion,seed,snrDb,sigma,truthMode,layer,truth,C,fit,calTime,fitTime,...
                        gLow,gHigh,deltaTrue,lowSeed,highSeed);
                catch ME
                    fitTime=toc(fitTic);row=pack_error(frameworkVersion,seed,snrDb,sigma,truthMode,layer,...
                        calTime,fitTime,lowSeed,highSeed,ME);
                end
                ir=ir+1;rows(ir)=row;Detail=struct2table(rows(1:ir));
                writetable(Detail,checkpoint);save(fullfile(outDir,'checkpoint.mat'),...
                    'Detail','ReferenceNoise','cfg0','caseModes','layers','-v7.3');
                fprintf('%s %s SNR=%g seed=%d -> %s, structure=%d, %.1fs\n',...
                    layer,truthMode,snrDb,seed,row.selected_mode,row.structure_correct,fitTime);
            end
        end
    end
end

Detail=struct2table(rows(1:ir));Summary=summarize_detail(Detail);
blind=Detail(Detail.validation_layer=="blind",:);Confusion=make_confusion(blind,caseModes);
writetable(Detail,fullfile(outDir,'detail.csv'));writetable(Summary,fullfile(outDir,'summary.csv'));
writetable(Confusion,fullfile(outDir,'blind_confusion_matrix.csv'));
save(fullfile(outDir,'hierarchical_structure_framework.mat'),'Detail','Summary','Confusion',...
    'ReferenceNoise','cfg0','caseModes','layers','frameworkVersion','-v7.3');
plot_framework_summary(Summary,outDir);plot_blind_confusion(Confusion,outDir);
disp(Summary);disp(Confusion);
end

function truth=case_truth(mode,fRot)
switch mode
    case "single_sync"
        truth=struct('f',10*fRot,'A',.25,'phi',pi/4);
    case "single_async"
        truth=struct('f',531,'A',.25,'phi',pi/4);
    case "dual_sync_sync"
        truth=struct('f',[10 26]*fRot,'A',[.25 .15],'phi',[pi/4 -pi/3]);
    case "dual_sync_async"
        truth=struct('f',[10*fRot 917.5],'A',[.25 .15],'phi',[pi/4 -pi/3]);
    otherwise
        error('Unknown truth mode: %s.',mode);
end
end

function sigma=reference_noise_std(model,delta,cfg,domain,snrList)
clean=simulate_highspeed_from_low_increment(model,delta,cfg,Inf,'snr_db',1);
q=cfg;q.A_true=[];q.f_true=[];q.phi_true=[];
stat=simulate_highspeed_from_low_increment(model,delta,q,Inf,'snr_db',1);
mv=map_highspeed_to_space(clean,cfg.alpha_k,cfg.R_tip,domain,.02);
ms=map_highspeed_to_space(stat,cfg.alpha_k,cfg.R_tip,domain,.02);
sigma=rms(mv.V_a-ms.V_a)./10.^(snrList/20);
end

function tf=is_completed(rows,n,seed,snr,mode,layer)
tf=false;
for i=1:n
    if rows(i).seed==seed&&rows(i).reference_snr_db==snr&&...
            string(rows(i).truth_mode)==mode&&string(rows(i).validation_layer)==layer
        tf=true;return
    end
end
end

function r=empty_row()
r=struct('framework_version',"",'seed',NaN,'reference_snr_db',NaN,'fixed_noise_std_V',NaN,...
    'low_seed',NaN,'high_seed',NaN,'truth_mode',"",'validation_layer',"",...
    'selected_mode',"",'structure_correct',false,'model_order_true',NaN,...
    'model_order_est',NaN,'f1_true_hz',NaN,'f2_true_hz',NaN,'f1_est_hz',NaN,...
    'f2_est_hz',NaN,'max_frequency_error_hz',NaN,'A1_true_mm',NaN,'A2_true_mm',NaN,...
    'A1_est_mm',NaN,'A2_est_mm',NaN,'max_amplitude_error_mm',NaN,...
    'g_low_est_mm',NaN,'g_high_est_mm',NaN,'delta_gap_est_mm',NaN,...
    'delta_gap_error_mm',NaN,'rmse_V',NaN,'selected_window_mm',NaN,...
    'calibration_time_s',NaN,'fit_time_s',NaN,'router_budget',"",...
    'router_bic_margin',NaN,'parameter_success',false,'error_id',"",'error_message',"");
end

function r=pack_row(version,seed,snr,sigma,truthMode,layer,truth,C,fit,calTime,fitTime,...
        gLow,gHigh,deltaTrue,lowSeed,highSeed)
r=empty_row();r.framework_version=version;r.seed=seed;r.reference_snr_db=snr;r.fixed_noise_std_V=sigma;
r.low_seed=lowSeed;r.high_seed=highSeed;r.truth_mode=truthMode;r.validation_layer=layer;
r.selected_mode=string(get_field(fit,'structure_selected',fit.mode));
r.structure_correct=r.selected_mode==truthMode;r.model_order_true=numel(truth.f);
r.model_order_est=fit.model_order;
[ft,At,fe,Ae]=match_components(truth.f,truth.A,fit.f_id,fit.A_id);
r.f1_true_hz=ft(1);r.f1_est_hz=fe(1);r.A1_true_mm=At(1);r.A1_est_mm=Ae(1);
if numel(ft)>1,r.f2_true_hz=ft(2);r.f2_est_hz=fe(2);r.A2_true_mm=At(2);r.A2_est_mm=Ae(2);end
r.max_frequency_error_hz=max(abs(fe-ft));r.max_amplitude_error_mm=max(abs(Ae-At));
r.g_low_est_mm=C.pathCal.g0;r.g_high_est_mm=fit.g_used;
r.delta_gap_est_mm=fit.g_used-C.pathCal.g0;r.delta_gap_error_mm=abs(r.delta_gap_est_mm-deltaTrue);
r.rmse_V=fit.rmse;r.selected_window_mm=get_field(C,'low_speed_selected_window_mm',NaN);
r.calibration_time_s=calTime;r.fit_time_s=fitTime;
r.router_budget=string(get_field(fit,'structure_budget',"known_structure"));
r.router_bic_margin=get_field(fit,'structure_bic_margin',NaN);
r.parameter_success=r.structure_correct&&r.max_amplitude_error_mm<=.01&&...
    r.delta_gap_error_mm<=.01&&r.max_frequency_error_hz<=2;
if abs(gLow-.5)>eps||abs(gHigh-.6)>eps,error('Unexpected gap truth.');end
end

function [ft,At,fe,Ae]=match_components(ft,At,fe,Ae)
ft=ft(:).';At=At(:).';fe=fe(:).';Ae=Ae(:).';
if numel(ft)~=numel(fe)
    fe=nan(size(ft));Ae=nan(size(At));return
end
if numel(ft)==2&&numel(fe)==2
    if sum(abs(fe([2 1])-ft))<sum(abs(fe-ft)),fe=fe([2 1]);Ae=Ae([2 1]);end
end
end

function r=pack_error(version,seed,snr,sigma,mode,layer,calTime,fitTime,lowSeed,highSeed,ME)
r=empty_row();r.framework_version=version;r.seed=seed;r.reference_snr_db=snr;r.fixed_noise_std_V=sigma;
r.truth_mode=mode;r.validation_layer=layer;r.low_seed=lowSeed;r.high_seed=highSeed;
r.calibration_time_s=calTime;r.fit_time_s=fitTime;r.error_id=string(ME.identifier);
r.error_message=string(ME.message);
end

function S=summarize_detail(D)
[G,layer,mode,snr]=findgroups(D.validation_layer,D.truth_mode,D.reference_snr_db);
S=table(layer,mode,snr,splitapply(@numel,D.seed,G),...
    splitapply(@mean,D.structure_correct,G),splitapply(@mean,D.parameter_success,G),...
    splitapply(@(x)median(x,'omitnan'),D.max_amplitude_error_mm,G),...
    splitapply(@(x)median(x,'omitnan'),D.delta_gap_error_mm,G),...
    splitapply(@(x)median(x,'omitnan'),D.fit_time_s,G),...
    'VariableNames',{'validation_layer','truth_mode','reference_snr_db','n',...
    'structure_correct_rate','parameter_success_rate','amplitude_error_median_mm',...
    'delta_gap_error_median_mm','fit_time_median_s'});
end

function C=make_confusion(D,modes)
M=zeros(numel(modes),numel(modes));
for i=1:numel(modes)
    for j=1:numel(modes)
        M(i,j)=sum(D.truth_mode==modes(i)&D.selected_mode==modes(j));
    end
end
C=array2table(M,'VariableNames',cellstr("est_"+modes));
C=addvars(C,modes(:),'Before',1,'NewVariableNames','truth_mode');
end

function plot_framework_summary(S,outDir)
if isempty(S),return;end
labels=categorical(S.validation_layer+" | "+S.truth_mode+" | "+...
    string(S.reference_snr_db)+" dB");
f=figure('Visible','off','Color','w','Position',[80 80 1500 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;bar(labels,[S.structure_correct_rate S.parameter_success_rate]);
ylim([0 1.05]);ylabel('Success rate');legend('Structure','Structure + parameters',...
    'Location','southoutside','Orientation','horizontal');grid on;title('Identification success');
nexttile;semilogy(labels,max(S.amplitude_error_median_mm,1e-7),'o-','LineWidth',1.3);hold on;
semilogy(labels,max(S.delta_gap_error_median_mm,1e-7),'s-','LineWidth',1.3);
yline(.01,'k--','0.01 mm');ylabel('Median absolute error (mm)');
legend('Amplitude','Clearance change','Tolerance','Location','best');grid on;title('Parameter errors');
nexttile([1 2]);bar(labels,S.fit_time_median_s);ylabel('Median time (s)');
grid on;title('Computational cost');
set(findall(f,'Type','axes'),'FontName','Times New Roman','FontSize',11);
exportgraphics(f,fullfile(outDir,'framework_summary.png'),'Resolution',220);close(f);
end

function plot_blind_confusion(C,outDir)
if isempty(C),return;end
M=C{:,2:end};if ~any(M,'all'),return;end
f=figure('Visible','off','Color','w','Position',[100 100 850 700]);
imagesc(M);axis image;colorbar;colormap(parula);
xticks(1:size(M,2));xticklabels(erase(string(C.Properties.VariableNames(2:end)),"est_"));
yticks(1:height(C));yticklabels(C.truth_mode);xlabel('Estimated structure');ylabel('True structure');
title('Blind structure confusion matrix');
for i=1:size(M,1)
    for j=1:size(M,2)
        text(j,i,string(M(i,j)),'HorizontalAlignment','center',...
            'FontWeight','bold','Color','w');
    end
end
set(gca,'FontName','Times New Roman','FontSize',11,'TickLabelInterpreter','none');
exportgraphics(f,fullfile(outDir,'blind_confusion_matrix.png'),'Resolution',220);close(f);
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
