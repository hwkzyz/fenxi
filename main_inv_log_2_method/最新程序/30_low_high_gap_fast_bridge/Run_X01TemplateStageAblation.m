function Stage=Run_X01TemplateStageAblation()
%RUN_X01TEMPLATESTAGEABLATION Identify the processing step that flips X01.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','x01_template_stage_ablation');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
cfg=make_cfg(make_inv_log_2_demo_case(ctx,.2).cfgCase);c=make_case();cfg.f_true=c.eo*cfg.RPM_high/60;cfg.A_true=c.A;cfg.phi_true=c.phi;
refModel=identity_model(lib,c.gLow,cfg);geom=simulate_highspeed_from_low_increment(refModel,.1,cfg,Inf,'fixed_std',1);
daq=make_daq(geom,lib,c,cfg);map=map_highspeed_to_space(daq,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
lowRaw=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),cfg,lib.domain,Inf,911001,'fixed_std');
formal=calibrate_inv_log_2_low_speed(lowRaw,lib,cfg);[Tb,xBin,count,stdv,sids,mapped]=raw_bins(lowRaw,cfg,lib);
span=formal.low_speed_selected_span_bins;
exact=repmat(eval_gap_template(lib,c.gLow,lib.xGrid),1,numel(cfg.alpha_k));
stages=cell(5,1);names=["exact reference";"bin only";"bin + SG";"bin + SG + PCHIP";"full bridge"];
stages{1}=refModel;
stages{2}=model_from_bins(Tb,xBin,lib,count,stdv,sids,mapped,'linear',1);
stages{3}=model_from_bins(Tb,xBin,lib,count,stdv,sids,mapped,'linear',span);
stages{4}=model_from_bins(Tb,xBin,lib,count,stdv,sids,mapped,'pchip',span);
stages{5}=formal.templateModel;
rows=repmat(row0(),numel(stages),1);
for i=1:numel(stages)
    m=stages{i};tf=fixed_pair(map,m,cfg,c.eo,c.gHigh);wf=fixed_pair(map,m,cfg,[10 18],c.gHigh);
    [vt,~]=predict(m,map,tf.theta,tf.eo);[vw,~]=predict(m,map,wf.theta,wf.eo);[vStage,~]=predict(m,map,true_theta(map,c,cfg),c.eo);
    [vRef,~]=predict(refModel,map,true_theta(map,c,cfg),c.eo);e=vStage-vRef;d=vw-vt;
    rows(i).stage=names(i);rows(i).template_rmse_V=sqrt(mean((m.templateLow-exact).^2,'all'));
    rows(i).true_rmse_V=tf.rmse;rows(i).wrong_rmse_V=wf.rmse;rows(i).delta_rmse_V=wf.rmse-tf.rmse;
    rows(i).winner=ternary(rows(i).delta_rmse_V>=0,"(10,12)","(10,18)");
    rows(i).eta_comp_V=abs(e'*d)/max(norm(d),eps);rows(i).rho_bias=abs(e'*d)/max(norm(e)*norm(d),eps);
    rows(i).g0_est_mm=m.pathCal.g0;
end
Stage=struct2table(rows);writetable(Stage,fullfile(outDir,'summary.csv'));save(fullfile(outDir,'stage_ablation.mat'),'Stage','formal','map','c','cfg','-v7.3');
write_report(outDir);disp(Stage);
end

function cfg=make_cfg(cfg)
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;cfg.route30ForwardModel="low_increment";
cfg.route30LowTemplateMethod="adaptive_sg";cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,'minCoverage',.90,'numFolds',5,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg.route30DualFrequencyRangeHz=[500 1300];cfg.route30GapHalfWidthMm=.20;cfg.dualSyncGapCount=11;cfg.funnelAllPairsFirstJointStep=true;cfg.funnelSecondJointPairCount=24;cfg.funnelFullRefineCount=1;cfg.dualSyncMaxIter=600;cfg.structuredComponentAmplitudeFloorMm=.05;cfg.returnDualSyncDiagnostics=true;
end

function c=make_case(),c=struct('eo',[10 12],'A',[.25 .10],'phi',[0 pi/2],'gLow',.5,'gHigh',.7);end

function m=identity_model(lib,g,cfg)
T=repmat(eval_gap_template(lib,g,lib.xGrid),1,numel(cfg.alpha_k));pc=struct('g0',g,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));m=build_path_template_model(lib,T,pc);
end

function [Tb,xBin,C,S,sids,mapped]=raw_bins(low,cfg,lib)
mapped=map_highspeed_to_space(low,cfg.alpha_k,cfg.R_tip,lib.domain,0);dx=.02;xBin=(lib.domain(1):dx:lib.domain(2)).';sids=unique(mapped.S_v(:),'stable');Tb=nan(numel(xBin),numel(sids));C=zeros(size(Tb));S=nan(size(Tb));
for is=1:numel(sids)
    q=mapped.S_v==sids(is);ib=round((mapped.x_v(q)-xBin(1))/dx)+1;v=mapped.V_a(q);ok=isfinite(v)&ib>=1&ib<=numel(xBin);ib=ib(ok);v=v(ok);
    for j=1:numel(xBin),z=v(ib==j);C(j,is)=numel(z);if numel(z)>=5,Tb(j,is)=mean(z);S(j,is)=std(z);end,end
    Tb(:,is)=fillmissing(Tb(:,is),'linear','EndValues','nearest');S(:,is)=fillmissing(S(:,is),'nearest');
end
end

function m=model_from_bins(Tb,xBin,lib,C,S,sids,mapped,interpMethod,span)
T=Tb;if span>=3,for j=1:size(T,2),T(:,j)=smoothdata(T(:,j),'sgolay',span);end,end
Tout=interp1(xBin,T,lib.xGrid,interpMethod,'extrap');low=struct('xGrid',lib.xGrid,'templateLow',mean(Tout,2),...
    'templateBySensor',Tout,'sensorIds',sids,'countLow',sum(interp1(xBin,C,lib.xGrid,'nearest','extrap'),2),...
    'stdLow',mean(interp1(xBin,S,lib.xGrid,'nearest','extrap'),2),'mapped',mapped);
opts=struct('xDomain',lib.domain,'fitMask',[true false false false],'fitGainOffset',false,'g0Init',.5,'sigmaFloor',1);
[pc,~]=calibrate_low_speed_path(lib,low,opts);pc.sensorIds=sids;m=build_path_template_model(lib,Tout,pc);
end

function D=make_daq(geom,lib,c,cfg)
t=geom.t(:);omega=cfg.RPM_high*2*pi/60;xHalf=.52*diff(lib.domain);u=c.A(1)*sin(2*pi*cfg.f_true(1)*t+c.phi(1))+c.A(2)*sin(2*pi*cfg.f_true(2)*t+c.phi(2));V=min(eval_gap_template(lib,c.gLow,lib.xGrid))*ones(size(t));Vs=V;
for rev=1:cfg.NumRevs_high,for s=1:numel(cfg.alpha_k),te=geom.T_opr_truth(rev)+cfg.alpha_k(s)/omega;idx=find(abs(t-te)<=xHalf/(omega*cfg.R_tip));x=omega*cfg.R_tip*(t(idx)-te);q=x-u(idx);ok=q>=lib.domain(1)&q<=lib.domain(2);if any(ok),V(idx(ok))=eval_gap_template(lib,c.gHigh,q(ok),s);end;ok=x>=lib.domain(1)&x<=lib.domain(2);if any(ok),Vs(idx(ok))=eval_gap_template(lib,c.gHigh,x(ok),s);end,end,end
D=geom;D.V_clean=V;D.V_static=Vs;D.V_cap=V;
end

function fit=fixed_pair(map,m,cfg,pair,g)
q=cfg;q.dualSyncCandidateEOPairs=pair;state=struct('gHat',g,'g_low_hat',g,'dx0',0);fit=run_dual_sync_voltage_vp_funnel(map,m,q,state).fit;
end

function z=true_theta(map,c,cfg)
off=angle(mean(exp(1i*(map.theta_v(:)-2*pi*(cfg.RPM_high/60)*map.t_v(:)))));p=c.phi-c.eo*off;z=[c.gHigh 0 c.A(1)*cos(p(1)) c.A(1)*sin(p(1)) c.A(2)*cos(p(2)) c.A(2)*sin(p(2))];
end

function [v,r]=predict(m,map,z,pair)
th=map.theta_v(:);u=z(3)*sin(pair(1)*th)+z(4)*cos(pair(1)*th)+z(5)*sin(pair(2)*th)+z(6)*cos(pair(2)*th);v=eval_gap_template(m,z(1),map.x_v(:)-z(2)-u,map.S_v(:));r=map.V_a(:)-v;v(~isfinite(v))=0;r(~isfinite(r))=0;
end

function r=row0(),r=struct('stage',"",'template_rmse_V',NaN,'true_rmse_V',NaN,'wrong_rmse_V',NaN,'delta_rmse_V',NaN,'winner',"",'eta_comp_V',NaN,'rho_bias',NaN,'g0_est_mm',NaN);end
function y=ternary(tf,a,b),if tf,y=a;else,y=b;end,end
function write_report(outDir)
fid=fopen(fullfile(outDir,'README.md'),'w');fprintf(fid,'# X01 Template-Stage Ablation\n\n`delta_rmse_V = RMSE_(10,18) - RMSE_(10,12)`: positive values select the true EO. `eta_comp_V` and `rho_bias` measure the stage-induced high-speed model error along the fitted EO discrimination direction.\n');fclose(fid);
end
