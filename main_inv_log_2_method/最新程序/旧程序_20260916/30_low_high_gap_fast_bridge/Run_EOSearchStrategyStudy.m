function [Detail,Summary]=Run_EOSearchStrategyStudy(numSeeds,snrList,smoothSpans)
%RUN_EOSEARCHSTRATEGYSTUDY Compare three unknown-EO search allocations.
% A: global top-N starts; B: one start for every EO pair; C: Top-M pairs
% followed by K nonlinear starts per pair. Low/high data and the final
% nonlinear residual are identical across strategies within each case.
if nargin<1||isempty(numSeeds),numSeeds=3;end
if nargin<2||isempty(snrList),snrList=[25 15 5];end
if nargin<3||isempty(smoothSpans),smoothSpans=[9 81];end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','eo_search_strategy_study');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg0=base.cfgCase;
cfg0.RPM_low=min(cfg0.RPM_high,300);cfg0.NumRevs_low=20;cfg0.NumRevs_high=8;
cfg0.route30ForwardModel="low_increment";cfg0.A_true=[.25 .15];
fRot=cfg0.RPM_high/60;eoTrue=[10 26];cfg0.f_true=eoTrue*fRot;cfg0.phi_true=[pi/4 -pi/3];
cfg0.route30DualFrequencyRangeHz=[min(cfg0.f_true) max(cfg0.f_true)];
cfg0.route30GapHalfWidthMm=.20;cfg0.dualSyncGapCount=11;cfg0.dualSyncKeepPerGap=80;
cfg0.dualSyncIteratedReplayCount=2100;cfg0.dualSyncRefineCount=50;
cfg0.dualSyncBalancedStartsPerPair=2;cfg0.dualSyncTopPairCount=40;
cfg0.dualSyncStartsPerTopPair=2;cfg0.structuredComponentAmplitudeFloorMm=.05;
cfg0.returnDualSyncDiagnostics=true;
strategyKeys=["global_budget","balanced_per_pair","coarse_top_m_multistart"];
strategyLabels=["A: global-50","B: stratified-2/pair","C: Top40 x 2"];
g0=.50;delta=.10;
truthT=repmat(eval_gap_template(lib,g0,lib.xGrid),1,numel(cfg0.alpha_k));
truthCal=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg0.alpha_k));
truthModel=build_path_template_model(lib,truthT,truthCal);
nMax=numSeeds*numel(snrList)*numel(smoothSpans)*numel(strategyKeys);
rows=repmat(empty_row(),nMax,1);ir=0;checkpointCsv=fullfile(outDir,'detail_checkpoint.csv');
if exist(checkpointCsv,'file')
    old=readtable(checkpointCsv,'TextType','string');
    keepOld=ismember(old.seed,1:numSeeds)&ismember(old.snr_db,snrList)&ismember(old.smooth_span,smoothSpans);
    keepStrategy=false(height(old),1);
    for ik=1:numel(strategyKeys)
        keepStrategy=keepStrategy|(old.strategy_key==strategyKeys(ik)&old.strategy_label==strategyLabels(ik));
    end
    keepOld=keepOld&keepStrategy;
    old=old(keepOld,:);oldRows=table2struct(old);nOld=min(numel(oldRows),nMax);
    rows(1:nOld)=oldRows(1:nOld);ir=nOld;
    fprintf('Resuming from %d completed strategy runs.\n',ir);
end
for isnr=1:numel(snrList)
    snrDb=snrList(isnr);
    for iseed=1:numSeeds
        cfgData=cfg0;
        lowData=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfgData,lib.domain,snrDb,941000+iseed);
        high=simulate_highspeed_from_low_increment(truthModel,delta,cfgData,snrDb,'snr_db',942000+iseed);
        highMap=map_highspeed_to_space(high,cfgData.alpha_k,cfgData.R_tip,lib.domain,.02);
        for isp=1:numel(smoothSpans)
            span=smoothSpans(isp);
            low=build_low_speed_templates_binned(lowData,cfgData.alpha_k,cfgData.R_tip,lib.domain,lib.xGrid,'mean',...
                struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',span,'minCoverage',.90));
            opts=struct('xDomain',lib.domain,'fitMask',[true false false false],...
                'fitGainOffset',false,'g0Init',g0,'sigmaFloor',1);
            [pc,cal]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=low.sensorIds;
            model=build_path_template_model(lib,low.templateBySensor,pc);
            state=estimate_highspeed_static_gap_raw(low.mapped,lib,cfgData);
            state.gHat=pc.g0;state.g_low_hat=pc.g0;state.dx0=0;
            for ist=1:numel(strategyKeys)
                if is_completed(rows,ir,iseed,snrDb,span,strategyKeys(ist),strategyLabels(ist)),continue;end
                cfg=cfgData;cfg.dualSyncSearchStrategy=strategyKeys(ist);
                runTic=tic;fit=run_dual_sync_voltage_vp(highMap,model,cfg,state);elapsed=toc(runTic);
                ir=ir+1;rows(ir)=pack_row(iseed,snrDb,span,strategyKeys(ist),strategyLabels(ist),...
                    low,truthT,lib.xGrid,pc,cal,fit,elapsed,g0,delta,cfg.A_true,eoTrue);
                Detail=struct2table(rows(1:ir));
                writetable(Detail,checkpointCsv);
                save(fullfile(outDir,'study_checkpoint.mat'),'Detail','cfg0','strategyKeys','strategyLabels',...
                    'snrList','smoothSpans','eoTrue','-v7.3');
                fprintf('SNR=%g seed=%d SG=%d %s -> EO=(%d,%d), true rank=%g, refined=%d, %.1fs\n',...
                    snrDb,iseed,span,strategyLabels(ist),fit.eo_id(1),fit.eo_id(2),...
                    rows(ir).true_pair_rank,rows(ir).true_pair_refined,elapsed);
            end
        end
    end
end
Detail=struct2table(rows(1:ir));Summary=summarize_detail(Detail,snrList,smoothSpans,strategyKeys,strategyLabels);
writetable(Detail,fullfile(outDir,'eo_search_strategy_detail.csv'));
writetable(Summary,fullfile(outDir,'eo_search_strategy_summary.csv'));
save(fullfile(outDir,'eo_search_strategy_study.mat'),'Detail','Summary','cfg0','strategyKeys',...
    'strategyLabels','snrList','smoothSpans','eoTrue','truthT','-v7.3');
plot_results(Summary,snrList,smoothSpans,strategyKeys,strategyLabels,outDir);
disp(Summary);
end

function tf=is_completed(rows,n,seed,snrDb,span,key,label)
tf=false;if n==0,return;end
for i=1:n
    if rows(i).seed==seed&&rows(i).snr_db==snrDb&&rows(i).smooth_span==span&&...
            string(rows(i).strategy_key)==key&&string(rows(i).strategy_label)==label
        tf=true;return
    end
end
end

function r=empty_row()
r=struct('seed',NaN,'snr_db',NaN,'smooth_span',NaN,'strategy_key',"",'strategy_label',"",...
    'eo1_est',NaN,'eo2_est',NaN,'eo_correct',false,'eo_max_error',NaN,...
    'A1_est_mm',NaN,'A2_est_mm',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,...
    'delta_gap_est_mm',NaN,'delta_gap_error_mm',NaN,'g0_error_mm',NaN,...
    'rmse_V',NaN,'true_pair_rank',NaN,'true_pair_refined',false,...
    'true_pair_refine_start_count',0,'true_pair_best_refined_rmse_V',NaN,...
    'template_rmse_V',NaN,'weighted_derivative_rmse_V_per_mm',NaN,...
    'calibration_rmse_V',NaN,'runtime_s',NaN,'solver_time_s',NaN,...
    'coarse_time_s',NaN,'replay_time_s',NaN,'refine_time_s',NaN,...
    'coarse_candidate_count',NaN,'iterated_replay_count',NaN,'nonlinear_refine_count',NaN);
end

function r=pack_row(seed,snrDb,span,key,label,low,Ttrue,x,pc,cal,fit,elapsed,g0,delta,Atrue,eoTrue)
r=empty_row();r.seed=seed;r.snr_db=snrDb;r.smooth_span=span;
r.strategy_key=key;r.strategy_label=label;eo=fit.eo_id(:).';A=fit.A_id(:).';
r.eo1_est=eo(1);r.eo2_est=eo(2);r.eo_correct=isequal(eo,eoTrue);
r.eo_max_error=max(abs(eo-eoTrue));r.A1_est_mm=A(1);r.A2_est_mm=A(2);
r.A1_error_mm=abs(A(1)-Atrue(1));r.A2_error_mm=abs(A(2)-Atrue(2));
r.delta_gap_est_mm=fit.g_used-pc.g0;r.delta_gap_error_mm=abs(r.delta_gap_est_mm-delta);
r.g0_error_mm=abs(pc.g0-g0);r.rmse_V=fit.rmse;r.calibration_rmse_V=cal.rmse_V;
d=median(diff(x));D=spatial_gradient(low.templateBySensor,d);Dt=spatial_gradient(Ttrue,d);w=Dt.^2;
r.template_rmse_V=sqrt(mean((low.templateBySensor-Ttrue).^2,'all'));
r.weighted_derivative_rmse_V_per_mm=sqrt(sum(w.*(D-Dt).^2,'all')/sum(w,'all'));
r.runtime_s=elapsed;r.solver_time_s=fit.solver_time_s;r.coarse_time_s=fit.coarse_time_s;
r.replay_time_s=fit.replay_time_s;r.refine_time_s=fit.refine_time_s;
r.coarse_candidate_count=fit.coarse_candidate_count;r.iterated_replay_count=fit.iterated_replay_count;
r.nonlinear_refine_count=fit.nonlinear_refine_count;
if isfield(fit,'diagnostic_pair_eo')
    q=find(all(fit.diagnostic_pair_eo==eoTrue,2),1);
    if ~isempty(q),r.true_pair_rank=fit.diagnostic_pair_rank(q);end
end
if isfield(fit,'diagnostic_refine_start_eo')
    r.true_pair_refine_start_count=sum(all(fit.diagnostic_refine_start_eo==eoTrue,2));
end
if isfield(fit,'diagnostic_refined_eo')
    q=all(fit.diagnostic_refined_eo==eoTrue,2);r.true_pair_refined=any(q);
    if any(q),r.true_pair_best_refined_rmse_V=min(fit.diagnostic_refined_rmse(q));end
end
end

function S=summarize_detail(D,snrList,spans,keys,labels)
n=numel(snrList)*numel(spans)*numel(keys);s=repmat(empty_summary(),n,1);ir=0;
for isp=1:numel(spans)
    for isnr=1:numel(snrList)
        for ik=1:numel(keys)
            q=D.smooth_span==spans(isp)&D.snr_db==snrList(isnr)&D.strategy_key==keys(ik);X=D(q,:);ir=ir+1;
            s(ir).smooth_span=spans(isp);s(ir).snr_db=snrList(isnr);s(ir).strategy_key=keys(ik);s(ir).strategy_label=labels(ik);
            s(ir).n=height(X);s(ir).eo_correct_rate=mean(X.eo_correct);s(ir).true_pair_refined_rate=mean(X.true_pair_refined);
            s(ir).true_pair_rank_median=median(X.true_pair_rank);s(ir).A1_error_mm_median=median(X.A1_error_mm);
            s(ir).A2_error_mm_median=median(X.A2_error_mm);s(ir).max_A_error_mm_median=median(max(X{:,{'A1_error_mm','A2_error_mm'}},[],2));
            s(ir).delta_gap_error_mm_median=median(X.delta_gap_error_mm);s(ir).true_pair_best_rmse_V_median=median(X.true_pair_best_refined_rmse_V,'omitnan');
            s(ir).final_rmse_V_median=median(X.rmse_V);s(ir).runtime_s_median=median(X.runtime_s);
            s(ir).nonlinear_refine_count_median=median(X.nonlinear_refine_count);s(ir).candidate_count_median=median(X.coarse_candidate_count);
        end
    end
end
S=struct2table(s);
end

function s=empty_summary()
s=struct('smooth_span',NaN,'snr_db',NaN,'strategy_key',"",'strategy_label',"",'n',0,...
    'eo_correct_rate',NaN,'true_pair_refined_rate',NaN,'true_pair_rank_median',NaN,...
    'A1_error_mm_median',NaN,'A2_error_mm_median',NaN,'max_A_error_mm_median',NaN,...
    'delta_gap_error_mm_median',NaN,'true_pair_best_rmse_V_median',NaN,...
    'final_rmse_V_median',NaN,'runtime_s_median',NaN,'nonlinear_refine_count_median',NaN,...
    'candidate_count_median',NaN);
end

function plot_results(S,snrList,spans,keys,labels,outDir)
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 13]);
tiledlayout(numel(spans),3,'TileSpacing','compact','Padding','compact');
cc=[.30 .48 .70;.85 .45 .30;.38 .68 .55];
barHandles=[];
for isp=1:numel(spans)
    span=spans(isp);Q=S(S.smooth_span==span,:);
    ax=nexttile;Y=metric_matrix(Q,snrList,keys,'eo_correct_rate');b=bar(1:numel(snrList),100*Y,'grouped');
    for k=1:numel(b),b(k).FaceColor=cc(k,:);end
    set(ax,'XTick',1:numel(snrList),'XTickLabel',string(snrList));ylim([0 105]);ylabel('Correct EO pair (%)');xlabel('Common SNR (dB)');title(sprintf('SG span = %d',span));
    text(ax,-.13,1.04,char('a'+(isp-1)*3),'Units','normalized','FontWeight','bold','FontSize',9);
    if isp==1,barHandles=b;end
    ax=nexttile;hold on;
    for k=1:numel(keys)
        q=Q.strategy_key==keys(k);[~,o]=ismember(snrList,Q.snr_db(q));v=Q.max_A_error_mm_median(q);v=v(o);
        plot((1:numel(snrList))+.045*(k-2),1000*v,'-o','Color',cc(k,:),'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor','w');
    end
    set(gca,'XTick',1:numel(snrList),'XTickLabel',string(snrList));xlabel('Common SNR (dB)');ylabel('Median max amplitude error (\mum)');
    text(ax,-.13,1.04,char('b'+(isp-1)*3),'Units','normalized','FontWeight','bold','FontSize',9);
    ax=nexttile;Y=metric_matrix(Q,snrList,keys,'runtime_s_median');b=bar(1:numel(snrList),Y,'grouped');
    for k=1:numel(b),b(k).FaceColor=cc(k,:);end
    set(ax,'XTick',1:numel(snrList),'XTickLabel',string(snrList));ylabel('Median solver time (s)');xlabel('Common SNR (dB)');
    text(ax,-.13,1.04,char('c'+(isp-1)*3),'Units','normalized','FontWeight','bold','FontSize',9);
end
set(findall(fig,'Type','axes'),'TickDir','in','Box','on','FontName','Times New Roman','FontSize',8,'LineWidth',.7);
lg=legend(barHandles,cellstr(labels),'Box','off','Orientation','horizontal');lg.Layout.Tile='south';
sgtitle('EO search allocation controls robustness independently of template smoothing','FontName','Times New Roman','FontSize',10,'FontWeight','bold');
exportgraphics(fig,fullfile(outDir,'eo_search_strategy_comparison.png'),'Resolution',300);
try,exportgraphics(fig,fullfile(outDir,'eo_search_strategy_comparison.pdf'),'ContentType','vector');catch,end
try,print(fig,fullfile(outDir,'eo_search_strategy_comparison.svg'),'-dsvg');catch,end
close(fig);
end

function Y=metric_matrix(Q,snrList,keys,field)
Y=nan(numel(snrList),numel(keys));
for i=1:numel(snrList),for k=1:numel(keys),q=Q.snr_db==snrList(i)&Q.strategy_key==keys(k);Y(i,k)=Q.(field)(q);end,end
end

function D=spatial_gradient(T,d)
if isvector(T),D=gradient(T,d);else,[~,D]=gradient(T,d,d);end
end
