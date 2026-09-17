function [Detail,Summary]=Run_LowTemplateDerivativeRegularizationStudy(snrDb,numSeeds,numLowRevs)
%RUN_LOWTEMPLATEDERIVATIVEREGULARIZATIONSTUDY SG span trade-off.
% Noise-only diagnostic; no OPR drift, gain drift, or outlier turns.
if nargin<1||isempty(snrDb),snrDb=5;end
if nargin<2||isempty(numSeeds),numSeeds=10;end
if nargin<3||isempty(numLowRevs),numLowRevs=20;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=numLowRevs;g0=.5;truth=repmat(eval_gap_template(lib,g0,lib.xGrid),1,numel(cfg.alpha_k));d=median(diff(lib.xGrid));Dt=spatial_gradient(truth,d);w=Dt.^2;spans=[5 7 9 13 17 21 31 41 51 61 81 101];n=numSeeds*numel(spans);seedCol=zeros(n,1);spanCol=zeros(n,1);vRmse=zeros(n,1);dRmse=zeros(n,1);corrD=zeros(n,1);gainD=zeros(n,1);ir=0;
for seed=1:numSeeds
    data=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,snrDb,971000+seed);
    for span=spans
        low=build_low_speed_templates_binned(data,cfg.alpha_k,cfg.R_tip,lib.domain,lib.xGrid,'mean',struct('gridSpacingMm',.02,'minBinCount',5,'smoothSpan',span,'minCoverage',.90));D=spatial_gradient(low.templateBySensor,d);ir=ir+1;seedCol(ir)=seed;spanCol(ir)=span;vRmse(ir)=sqrt(mean((low.templateBySensor-truth).^2,'all'));dRmse(ir)=sqrt(sum(w.*(D-Dt).^2,'all')/sum(w,'all'));corrD(ir)=sum(D.*Dt,'all')/sqrt(sum(D.^2,'all')*sum(Dt.^2,'all'));gainD(ir)=sum(D.*Dt,'all')/sum(Dt.^2,'all');
    end
end
Detail=table(seedCol,spanCol,vRmse,dRmse,corrD,gainD,'VariableNames',{'seed','smooth_span','template_rmse_V','weighted_derivative_rmse_V_per_mm','derivative_correlation','derivative_gain'});Summary=table(spans(:),'VariableNames',{'smooth_span'});vars=Detail.Properties.VariableNames(3:end);for iv=1:numel(vars),v=vars{iv};med=zeros(numel(spans),1);iq=zeros(numel(spans),1);for i=1:numel(spans),q=Detail.(v)(Detail.smooth_span==spans(i));med(i)=median(q);iq(i)=iqr(q);end;Summary.([v,'_median'])=med;Summary.([v,'_iqr'])=iq;end
outDir=fullfile(root,'output','low_template_derivative_regularization');if ~exist(outDir,'dir'),mkdir(outDir);end;tag=sprintf('%gdb_%dseeds_%drevs',snrDb,numSeeds,numLowRevs);writetable(Detail,fullfile(outDir,['detail_',tag,'.csv']));writetable(Summary,fullfile(outDir,['summary_',tag,'.csv']));save(fullfile(outDir,['study_',tag,'.mat']),'Detail','Summary','spans','truth','lib');fig=figure('Color','w','Units','centimeters','Position',[2 2 17 7]);tiledlayout(1,3,'TileSpacing','compact','Padding','compact');fields={'template_rmse_V_median','weighted_derivative_rmse_V_per_mm_median','derivative_correlation_median'};iqfields={'template_rmse_V_iqr','weighted_derivative_rmse_V_per_mm_iqr','derivative_correlation_iqr'};scales=[1000 1 1];ylabs={'Template RMSE (mV)','Weighted derivative RMSE (V/mm)','Derivative correlation'};for p=1:3,nexttile;errorbar(Summary.smooth_span,scales(p)*Summary.(fields{p}),.5*scales(p)*Summary.(iqfields{p}),'-o','LineWidth',1.1,'MarkerSize',4);xlabel('SG smoothing span (bins)');ylabel(ylabs{p});box on;end;set(findall(fig,'Type','axes'),'TickDir','in','FontName','Times New Roman','FontSize',8);exportgraphics(fig,fullfile(outDir,['span_tradeoff_',tag,'.png']),'Resolution',240);try,exportgraphics(fig,fullfile(outDir,['span_tradeoff_',tag,'.pdf']),'ContentType','vector');catch,end;disp(Summary);
end
function D=spatial_gradient(T,d)
if isvector(T),D=gradient(T,d);else,[~,D]=gradient(T,d,d);end
end
