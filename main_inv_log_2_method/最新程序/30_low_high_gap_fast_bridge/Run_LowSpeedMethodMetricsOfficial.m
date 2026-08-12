function T=Run_LowSpeedMethodMetricsOfficial(snrList,lowSeeds)
%RUN_LOWSPEEDMETHODMETRICSOFFICIAL Low-speed-only metrics on the formal path.
% The record, support ratio, spatial grid, and reference curve match the
% official dual-sync generalization driver. High-speed data are not used.
if nargin<1||isempty(snrList),snrList=[25 15 5];end
if nargin<2||isempty(lowSeeds),lowSeeds=1:10;end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;
cfg.route30ForwardModel="low_increment";cfg.route30SupportAware=true;
cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm',...
    [.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
g0=.50;truth=eval_gap_template(lib,g0,lib.xGrid);truth=repmat(truth,1,numel(cfg.alpha_k));
dx=median(diff(lib.xGrid));Dtruth=gradmat(truth,dx);
methods=["M0_bin_mean";"M1_adaptive_sg";"M2_local_quadratic";"M3_spline"];
rows=repmat(empty_row(),numel(snrList)*numel(lowSeeds)*numel(methods),1);ir=0;
for snr=snrList(:).'
    for iseed=1:numel(lowSeeds)
        lowSeed=991001+round(1000*g0)+lowSeeds(iseed);
        low=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,Inf,...
            lowSeed,'fixed_std');
        sig=rms(low.V_clean-mean(low.V_clean));sigma=sig/10^(snr/20);
        rng(lowSeed+10000+snr,'twister');low.V_cap=low.V_clean+sigma*randn(size(low.V_clean));
        opt=cfg.route30LowTemplateOptions;opt.mapWindowRatio=low.spatial_window_ratio;
        for im=1:numel(methods)
        switch im
            case 1
                opt.smoothSpan=1;L=build_low_speed_templates_binned(low,cfg.alpha_k,cfg.R_tip,...
                    lib.domain,lib.xGrid,'mean',opt);
            case 2
                L=build_low_speed_templates_adaptive(low,cfg.alpha_k,cfg.R_tip,lib.domain,...
                    lib.xGrid,'cv_sg',opt);
            case 3
                L=build_low_speed_templates_local_quadratic(low,cfg.alpha_k,cfg.R_tip,...
                    lib.domain,lib.xGrid,opt);
            case 4
                L=build_low_speed_templates_adaptive(low,cfg.alpha_k,cfg.R_tip,lib.domain,...
                    lib.xGrid,'cv_spline',opt);
        end
        if ~isfield(L,'derivativeBySensor')
            L.derivativeBySensor=gradmat(L.templateBySensor,dx);
        end
            ir=ir+1;rows(ir)=pack_row(snr,lowSeeds(iseed),methods(im),L,truth,Dtruth,lib.xGrid);
        end
    end
end
T=struct2table(rows);outDir=fullfile(root,'output','official_low_speed_method_metrics');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(T,fullfile(outDir,'detail.csv'));S=groupsummary(T,'method',{'mean','median'},...
    {'template_rmse_V','derivative_rmse_V_per_mm','derivative_correlation',...
    'peak_bias_pct','width_bias_pct','area_bias_pct','slope_bias_pct'});
writetable(S,fullfile(outDir,'summary.csv'));disp(T);disp(S);
end

function r=empty_row()
r=struct('snr_db',NaN,'low_seed',NaN,'method',"",'template_rmse_V',NaN,...
    'derivative_rmse_V_per_mm',NaN,'derivative_correlation',NaN,...
    'peak_bias_pct',NaN,'width_bias_pct',NaN,'area_bias_pct',NaN,'slope_bias_pct',NaN);
end
function r=pack_row(snr,seed,m,L,T,D,x)
dx=median(diff(x));r=empty_row();r.snr_db=snr;r.method=m;
r.low_seed=seed;
r.template_rmse_V=sqrt(mean((L.templateBySensor-T).^2,'all'));
Q=L.derivativeBySensor;r.derivative_rmse_V_per_mm=sqrt(mean((Q-D).^2,'all'));
r.derivative_correlation=sum(Q.*D,'all')/sqrt(sum(Q.^2,'all')*sum(D.^2,'all'));
[p0,w0,a0,s0]=shape(T,D,x);[p1,w1,a1,s1]=shape(L.templateBySensor,Q,x);
r.peak_bias_pct=100*mean((p1-p0)./max(abs(p0),eps));
r.width_bias_pct=100*mean((w1-w0)./max(abs(w0),eps));
r.area_bias_pct=100*mean((a1-a0)./max(abs(a0),eps));
r.slope_bias_pct=100*mean((s1-s0)./max(abs(s0),eps));
end
function [pk,wd,ar,sl]=shape(T,D,x)
n=size(T,1);edge=max(3,round(.1*n));dx=median(diff(x));
pk=zeros(1,size(T,2));wd=pk;ar=pk;sl=pk;
for i=1:size(T,2)
    b=median([T(1:edge,i);T(end-edge+1:end,i)]);q=abs(T(:,i)-b);
    pk(i)=max(q);wd(i)=nnz(q>=.5*pk(i))*dx;ar(i)=trapz(x,q);sl(i)=max(abs(D(:,i)));
end
end
function D=gradmat(T,dx)
D=zeros(size(T));for i=1:size(T,2),D(:,i)=gradient(T(:,i),dx);end
end
