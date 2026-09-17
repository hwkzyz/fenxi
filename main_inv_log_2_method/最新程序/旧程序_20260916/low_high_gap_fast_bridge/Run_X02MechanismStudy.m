function Out=Run_X02MechanismStudy()
%RUN_X02MECHANISMSTUDY Clean-oracle and weak-EO competition diagnostics.
% This is a post-freeze analysis; it does not alter Funnel budgets.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','x02_mechanism_study');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
cfg=make_cfg(make_inv_log_2_demo_case(ctx,.2).cfgCase);c=struct('eo',[15 18],'A',[.15 .10],...
    'phi',[0 pi],'gLow',.5,'gHigh',.4);fRot=cfg.RPM_high/60;cfg.f_true=c.eo*fRot;cfg.A_true=c.A;cfg.phi_true=c.phi;
low=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),cfg,lib.domain,Inf,983001,'fixed_std');
cal=calibrate_inv_log_2_low_speed(low,lib,cfg);
truth=make_truth_model(lib,cfg,c.gLow);geom=simulate_highspeed_from_low_increment(truth,c.gHigh-c.gLow,cfg,Inf,'fixed_std',1);
daq=make_daq(geom,lib,c,cfg);map=map_highspeed_to_space(daq,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
lowProxy=struct('x_v',reshape(cal.low_support_observation.per_sensor_observed_mm.',[],1),...
    'S_v',repelem(cal.low_support_observation.sensor_ids(:),2));contract=derive_support_aware_domain(map,lowProxy,cal.templateModel,cfg.route30PhysicalDomain,cfg);
[map,contract]=restrict_map_to_support_contract(map,contract);cfg.supportContract=contract;
fit=run_dual_sync_voltage_vp_funnel(map,cal.templateModel,cfg,cal.lowState);fit.support_contract=contract;
[~,ord]=sort(fit.diagnostic_refined_rmse);
topEo=fit.diagnostic_refined_eo(ord(1:min(3,end)),:);topJ=fit.diagnostic_refined_rmse(ord(1:min(3,end)));
oracle=repmat(row0(),size(topEo,1)+1,1);pairs=[c.eo;topEo];
for i=1:size(pairs,1)
    q=cfg;q.dualSyncCandidateEOPairs=pairs(i,:);q.funnelFullRefineCount=1;
    q.supportContract=contract;r=run_dual_sync_voltage_vp_funnel(map,cal.templateModel,q,cal.lowState);oracle(i).eo1=pairs(i,1);oracle(i).eo2=pairs(i,2);
    oracle(i).rmse_V=r.rmse;oracle(i).is_true=all(pairs(i,:)==c.eo);oracle(i).VFit=r.VFit;
end
oracleTable=removevars(struct2table(oracle),{'VFit'});
oracleTable.delta_rmse_vs_true_V=oracleTable.rmse_V-oracleTable.rmse_V(1);
oracleTable.delta_J_vs_true_V2=oracleTable.rmse_V.^2-oracleTable.rmse_V(1).^2;
competitors=unique(topEo,'rows','stable');subspace=repmat(subspace_row(),size(competitors,1),1);
for i=1:size(competitors,1)
    subspace(i)=subspace_metric(fit,cfg,lib,map,competitors(i,:),c.eo(1));
end
snr=[25 15 5].';sigHigh=rms(daq.V_clean-mean(daq.V_clean));ic=find(~[oracle.is_true],1);if isempty(ic),ic=1;end;D=abs(oracle(1).VFit-oracle(ic).VFit);
metrics=table(snr,repmat(sqrt(mean(D.^2)),3,1),sqrt(mean(D.^2))./(sigHigh./10.^(snr/20)),...
    'VariableNames',{'nominal_snr_db','D_EO_V','Gamma_EO'});
Out=struct('oracle',oracleTable,'subspace',struct2table(subspace),'metrics',metrics,...
    'top_eo',topEo,'top_rmse_V',topJ,'support_contract',fit.support_contract);
writetable(oracleTable,fullfile(outDir,'oracle.csv'));writetable(Out.subspace,fullfile(outDir,'subspace.csv'));writetable(metrics,fullfile(outDir,'selection_distance.csv'));
save(fullfile(outDir,'study.mat'),'Out','fit','cal','map','-v7.3');disp(oracleTable);disp(Out.subspace);disp(metrics);
end
function cfg=make_cfg(cfg)
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;cfg.route30ForwardModel="low_increment";
cfg.route30LowTemplateMethod="adaptive_sg";cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,'minCoverage',.90,'numFolds',5,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg.route30DualFrequencyRangeHz=[500 1300];cfg.route30GapHalfWidthMm=.20;cfg.dualSyncGapCount=11;cfg.funnelAllPairsFirstJointStep=true;cfg.funnelSecondJointPairCount=24;cfg.funnelFullRefineCount=3;cfg.dualSyncMaxIter=600;cfg.structuredComponentAmplitudeFloorMm=.05;cfg.returnDualSyncDiagnostics=true;
cfg.route30SupportAware=true;cfg.route30PhysicalDomain=struct('dx_bounds_mm',[-.03 .03],'amplitude_max_mm',[.25 .25],'gap_bounds_mm',[.40 .70],'source',"predeclared paper simulation envelope");
end
function model=make_truth_model(lib,cfg,g)
T=repmat(eval_gap_template(lib,g,lib.xGrid),1,numel(cfg.alpha_k));pc=struct('g0',g,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));model=build_path_template_model(lib,T,pc);
end
function D=make_daq(G,lib,c,cfg)
t=G.t(:);om=cfg.RPM_high*2*pi/60;h=.52*diff(lib.domain);u=zeros(size(t));for j=1:2,u=u+c.A(j)*sin(2*pi*cfg.f_true(j)*t+c.phi(j));end
V=min(eval_gap_template(lib,c.gLow,lib.xGrid))*ones(size(t));for r=1:cfg.NumRevs_high,for s=1:numel(cfg.alpha_k),te=G.T_opr_truth(r)+cfg.alpha_k(s)/om;idx=find(abs(t-te)<=h/(om*cfg.R_tip));x=om*cfg.R_tip*(t(idx)-te);q=x-u(idx);ok=q>=lib.domain(1)&q<=lib.domain(2);if any(ok),V(idx(ok))=eval_gap_template(lib,c.gHigh,q(ok),s);end,end,end;D=G;D.V_clean=V;D.V_cap=V;
end
function r=subspace_metric(fit,~,lib,map,eo,common)
theta=map.theta_v(:);x=map.x_v(:);fx=eval_gap_derivative(lib,fit.g_used,x);h=.002;fg=(eval_gap_template(lib,fit.g_used+h,x)-eval_gap_template(lib,fit.g_used-h,x))/(2*h);
N=[-fx,fg,-fx.*sin(common*theta),-fx.*cos(common*theta)];[QN,~]=qr(N,0);
Btrue=unique_basis([15 18],common,fx,theta);B=unique_basis(eo,common,fx,theta);Btrue=Btrue-QN*(QN'*Btrue);B=B-QN*(QN'*B);[QT,~]=qr(Btrue,0);[QB,~]=qr(B,0);s=svd(QT'*QB);r=subspace_row();r.eo1=eo(1);r.eo2=eo(2);
if isempty(QB),r.rho_max=NaN;r.min_angle_deg=NaN;else,r.rho_max=max(s);r.min_angle_deg=acosd(min(1,max(s)));end
end
function B=unique_basis(pair,common,fx,theta)
u=pair(pair~=common);if isempty(u),B=zeros(numel(theta),0);return;end
B=[-fx.*sin(u(1)*theta),-fx.*cos(u(1)*theta)];
end
function r=row0(),r=struct('eo1',NaN,'eo2',NaN,'rmse_V',NaN,'is_true',false,'VFit',[]);end
function r=subspace_row(),r=struct('eo1',NaN,'eo2',NaN,'rho_max',NaN,'min_angle_deg',NaN);end
