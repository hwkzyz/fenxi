function T=Run_R2_CurvatureAudit()
% Compare empirical profile curvature with the nuisance-profiled Jacobian curvature.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');addpath(fullfile(mainDir,'identifiability_mechanism_analysis','local_func'),'-begin');
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
b=make_inv_log_2_demo_case(ctx,.2);cfg=b.cfgCase;cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;cfg.route30ForwardModel="low_increment";cfg.route30LowTemplateMethod="adaptive_sg";cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,'minCoverage',.90,'numFolds',5,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22]);
g0=.5;gt=.575;eo=10;A=.25;low=simulate_low_speed_template(@(z)eval_gap_template(lib,g0,z),cfg,lib.domain,Inf,99601,'fixed_std');C=calibrate_inv_log_2_low_speed(low,lib,cfg);d=simulate_rotating_waveform_from_template(@(z)eval_gap_template(lib,gt,z),cfg.RPM_high,cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,lib.domain,A,eo*cfg.RPM_high/60,pi/4,'noise_ratio');d.V_cap=d.V_clean;m=map_highspeed_to_space(d,cfg.alpha_k,cfg.R_tip,lib.domain,.02,.30);model=C.templateModel;
fun=@(z)res(z,gt,eo,m,model);[z,~]=solve_lsq_bounded(fun,[0 A/sqrt(2) A/sqrt(2)],[-.5 -1.5 -1.5],[.5 1.5 1.5],160,1e-11,1e-11);r0=fun(z);n=numel(r0);eg=1e-5;ez=1e-5;
Jg=(res(z,gt+eg,eo,m,model)-res(z,gt-eg,eo,m,model))/(2*eg);Jz=zeros(n,3);for j=1:3,zp=z;zm=z;zp(j)=zp(j)+ez;zm(j)=zm(j)-ez;Jz(:,j)=(res(zp,gt,eo,m,model)-res(zm,gt,eo,m,model))/(2*ez);end
[Q,~]=qr(Jz,0);qg=Jg-Q*(Q'*Jg);kJ=2*(qg'*qg)/n;
P=readtable(fullfile(root,'output','r2_heldout_profile.csv'));q=P.case_id==2&P.profile_type=="adaptive_profile";P=P(q,:);[~,i0]=min(abs(P.g_candidate-gt));h=P.g_candidate(i0+1)-P.g_candidate(i0);phi=P.rmse.^2;kEmp=(phi(i0+1)-2*phi(i0)+phi(i0-1))/h^2;
T=table(kEmp,kJ,kEmp/kJ,norm(Q' * Jg)^2/max(norm(Jg)^2,eps),n,'VariableNames',{'empirical_curvature','jacobian_profiled_curvature','curvature_ratio','gap_energy_absorbed_fraction','n_samples'});writetable(T,fullfile(root,'output','r2_curvature_audit.csv'));disp(T);
end
function r=res(z,g,eo,m,model),x=m.x_v(:);V=m.V_a(:);th=m.theta_v(:);r=V-eval_gap_template(model,g,x-z(1)-z(2)*sin(eo*th)-z(3)*cos(eo*th));r(~isfinite(r))=10*max(std(V),1e-3);end
