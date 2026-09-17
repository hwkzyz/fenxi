function Result=Run_BoundaryFailureOracleAudit()
%RUN_BOUNDARYFAILUREORACLEAUDIT Separate candidate failure from objective ambiguity.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
base=make_inv_log_2_demo_case(ctx,.2);cfg=base.cfgCase;
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=8;cfg.route30ForwardModel="absolute";
cfg.route30UseAllTrustedSamples=true;cfg.route30HierarchicalOrderSelection=true;
cfg.route30EnableRhoFallback=false;cfg.f1Grid=300:5:1500;cfg.f2Grid=cfg.f1Grid;
gLow=.8;gHigh=.5;phi=[pi/4 -pi/3];
cases={...
 struct('k',3,'name',"sep100_snr10",'snr',10,'f',[700 800],'A',[.25 .15]),...
 struct('k',5,'name',"sep50_snr5",'snr',5,'f',[700 750],'A',[.25 .15]),...
 struct('k',7,'name',"weak_snr5",'snr',5,'f',[700 1200],'A',[.25 .02]),...
 struct('k',8,'name',"weak_snr0",'snr',0,'f',[700 1200],'A',[.25 .02]),...
 struct('k',10,'name',"single_snr0",'snr',0,'f',531,'A',.25)};
rows=repmat(row0(),numel(cases),1);
for i=1:numel(cases)
    c=cases{i};cfg.snrDb=c.snr;
    low=simulate_low_speed_template(@(z)eval_gap_template(lib,gLow,z),...
        cfg,lib.domain,c.snr,203000+c.k);
    d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gHigh,z),...
        cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,...
        c.A,c.f,phi(1:numel(c.f)),'noise_ratio');
    sig=rms(d.V_clean-min(d.V_clean));rng(203100+c.k,'twister');
    d.V_cap=d.V_clean+sig/10^(c.snr/20)*randn(size(d.V_clean));
    map=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
    fit=run_inv_log_2_low_high_main(low,map,lib,cfg);
    rows(i).name=c.name;rows(i).selected_order=fit.model_order;rows(i).selected_rmse=fit.rmse;
    rows(i).single_rmse=fit.singleFit.rmse;rows(i).selected_f=join(string(fit.f_id),'+');
    if numel(c.f)==2
        [~,~,gid]=unique([map.rev_v(:),map.S_v(:)],'rows');
        p0=[c.A(1),phi(1),c.f(1),c.A(2),phi(2),c.f(2)];
        oracle=refine_projected_joint(map.t_v(:),map.V_a(:),map.x_v(:),lib,cfg,...
            gid,max(gid),gHigh,0,p0,0);
        Vo=eval_gap_template(lib,oracle.g,map.x_v(:)-oracle.dx-fit_u(oracle.p,map.t_v(:)));
        rows(i).oracle_rmse=sqrt(mean((map.V_a(:)-Vo).^2));
        rows(i).oracle_f=join(string(oracle.f),'+');
        rows(i).selected_minus_oracle_rmse=fit.rmse-rows(i).oracle_rmse;
    else
        u=c.A*sin(2*pi*c.f*map.t_v(:)+phi(1));
        Vt=eval_gap_template(lib,gHigh,map.x_v(:)-u);
        rows(i).oracle_rmse=sqrt(mean((map.V_a(:)-Vt).^2));rows(i).oracle_f=string(c.f);
        rows(i).selected_minus_oracle_rmse=fit.rmse-rows(i).oracle_rmse;
    end
    fprintf('%s: selected %s %.5g V, oracle %s %.5g V\n',c.name,...
        rows(i).selected_f,fit.rmse,rows(i).oracle_f,rows(i).oracle_rmse);
end
Audit=struct2table(rows);out=fullfile(root,'output','boundary_failure_oracle');
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'boundary_failure_oracle.csv'));save(fullfile(out,'boundary_failure_oracle.mat'),'Audit');disp(Audit);
Result=struct('Audit',Audit,'outputDir',out);
end

function r=row0()
r=struct('name',"",'selected_order',NaN,'selected_f',"",'selected_rmse',NaN,...
    'single_rmse',NaN,'oracle_f',"",'oracle_rmse',NaN,...
    'selected_minus_oracle_rmse',NaN);
end
