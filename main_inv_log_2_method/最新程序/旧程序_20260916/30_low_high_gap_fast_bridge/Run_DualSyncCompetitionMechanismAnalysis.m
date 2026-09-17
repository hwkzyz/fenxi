function [Detail,Summary,Subspace]=Run_DualSyncCompetitionMechanismAnalysis(lowSeeds,highSeeds,snrLevels)
%RUN_DUALSYNCCOMPETITIONMECHANISMANALYSIS Diagnose final EO competition.
% The Funnel V2 configuration remains fixed.  This study separates low-speed
% template noise from high-speed acquisition noise and compares the true X01
% pair (10,12) with its dominant observed competitor (10,18).
if nargin<1||isempty(lowSeeds),lowSeeds=1:5;end
if nargin<2||isempty(highSeeds),highSeeds=1:10;end
if nargin<3||isempty(snrLevels),snrLevels=[25 15 5];end

root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','dual_sync_competition_mechanism_v2');
if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
cfg0=make_cfg(make_inv_log_2_demo_case(ctx,.2).cfgCase);
conditions=make_conditions();
baseline=conditions(conditions.condition_id=="C01_baseline",:);
cfgBase=apply_condition(cfg0,baseline,cfg0.RPM_high/60);
truthBase=make_truth_model(lib,cfg0,baseline.g_low_mm);
geometry=simulate_highspeed_from_low_increment(truthBase,...
    baseline.g_high_mm-baseline.g_low_mm,cfgBase,Inf,'fixed_std',1);

% One clean/noise-source control for all pressure cases, then an X01 crossed
% low-template/high-record experiment.  In the two one-sided cells, a seed
% belonging to the clean source is intentionally ignored.
ids=["C02_eo_10_12";"X01_corner_near";"X02_corner_lowamp"];
pressure=conditions(ismember(conditions.condition_id,ids),:);
modes=["noise_free";"low_clean_high_noisy";"low_noisy_high_clean";"both_noisy"];
Detail=table();
for ic=1:height(pressure)
    for snr=snrLevels(:).'
        for im=1:numel(modes)
            Detail=[Detail;struct2table(run_cell(pressure(ic,:),cfg0,lib,geometry,snr,modes(im),1,1,"control",false))]; %#ok<AGROW>
        end
    end
end

x01=pressure(pressure.condition_id=="X01_corner_near",:);
for snr=snrLevels(:).'
    for il=1:numel(lowSeeds)
        for ih=1:numel(highSeeds)
            Detail=[Detail;struct2table(run_cell(x01,cfg0,lib,geometry,snr,"both_noisy",lowSeeds(il),highSeeds(ih),"crossed",false))]; %#ok<AGROW>
        end
    end
end

% The one-sided crossed cells make the attribution identifiable without
% needlessly repeating a source that is clean.
for snr=snrLevels(:).'
    for ih=1:numel(highSeeds)
        Detail=[Detail;struct2table(run_cell(x01,cfg0,lib,geometry,snr,"low_clean_high_noisy",1,highSeeds(ih),"crossed",false))]; %#ok<AGROW>
    end
    for il=1:numel(lowSeeds)
        Detail=[Detail;struct2table(run_cell(x01,cfg0,lib,geometry,snr,"low_noisy_high_clean",lowSeeds(il),1,"crossed",false))]; %#ok<AGROW>
    end
end

% Confirm the basin comparison itself with a deliberately expensive,
% multi-start optimization in clean and fully noisy controls at 25 dB.
for im=[1 4]
    Detail=[Detail;struct2table(run_cell(x01,cfg0,lib,geometry,25,modes(im),1,1,"oracle_multistart",true))]; %#ok<AGROW>
end

Summary=summarize_detail(Detail);
Subspace=build_subspace_audit(Detail);
writetable(Detail,fullfile(outDir,'detail.csv'));
writetable(Summary,fullfile(outDir,'summary.csv'));
writetable(Subspace,fullfile(outDir,'subspace.csv'));
save(fullfile(outDir,'competition_mechanism.mat'),'Detail','Summary','Subspace',...
    'lowSeeds','highSeeds','snrLevels','-v7.3');
write_report(Summary,Subspace,outDir);
disp(Summary);
end

function row=run_cell(c,cfg0,lib,geometry,snr,mode,lowIndex,highIndex,scope,highBudget)
cfg=apply_condition(cfg0,c,cfg0.RPM_high/60);
daq=make_condition_daq(geometry,lib,c,cfg,lib.domain);
[cal,lowMeta]=make_calibration(c,cfg0,lib,snr,mode,lowIndex);
[map,highMeta]=make_high_map(daq,cfg0,lib,snr,mode,highIndex);
fit=run_inv_log_2_high_with_calibration(cal,map,cfg);
truePair=[c.eo1_true c.eo2_true];wrongPair=[10 18];
[trueFit,trueMeta]=oracle_pair(map,cal.templateModel,cfg,fit,truePair,highBudget);
[wrongFit,wrongMeta]=oracle_pair(map,cal.templateModel,cfg,fit,wrongPair,highBudget);
subspace=subspace_metrics(map,cal.templateModel,trueFit,wrongFit);
row=pack_row(c,snr,mode,scope,lowIndex,highIndex,lowMeta,highMeta,fit,trueFit,wrongFit,trueMeta,wrongMeta,subspace);
end

function cfg=make_cfg(cfg)
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;
cfg.route30ForwardModel="low_increment";cfg.route30LowTemplateMethod="adaptive_sg";
cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,...
    'minCoverage',.90,'numFolds',5,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg.route30FrequencyStructureMode="dual_sync_sync";cfg.route30DualFrequencyRangeHz=[500 1300];
cfg.route30GapHalfWidthMm=.20;cfg.dualSyncGapCount=11;cfg.funnelFullRefineCount=3;
cfg.funnelAllPairsFirstJointStep=true;cfg.funnelSecondJointPairCount=24;
cfg.structuredComponentAmplitudeFloorMm=.05;cfg.returnDualSyncDiagnostics=true;
cfg.dualSyncSearchStrategy="funnel_v1";cfg.dualSyncMaxIter=180;
end

function [cal,meta]=make_calibration(c,cfg,lib,snr,mode,lowIndex)
seed=810001+10000*c.condition_seed_id+100*round(snr)+lowIndex;
low=simulate_low_speed_template(@(z)eval_gap_template(lib,c.g_low_mm,z),cfg,lib.domain,Inf,seed,'fixed_std');
isNoisy=(mode=="low_noisy_high_clean"||mode=="both_noisy");
if isNoisy
    signal=rms(low.V_clean-mean(low.V_clean));sigma=signal/10^(snr/20);rng(seed,'twister');
    low.V_cap=low.V_clean+sigma*randn(size(low.V_clean));
else
    signal=rms(low.V_clean-mean(low.V_clean));sigma=0;low.V_cap=low.V_clean;
end
cal=calibrate_inv_log_2_low_speed(low,lib,cfg);
meta=struct('seed',seed,'signal_rms_V',signal,'noise_std_V',sigma,'g0_est_mm',cal.pathCal.g0,...
    'selected_window_mm',get_field(cal,'low_speed_selected_window_mm',NaN));
end

function [map,meta]=make_high_map(daq,cfg,lib,snr,mode,highIndex)
seed=820001+10000*round(snr)+highIndex;
isNoisy=(mode=="low_clean_high_noisy"||mode=="both_noisy");
signal=rms(daq.V_clean-mean(daq.V_clean));
noisy=daq;
if isNoisy
    sigma=signal/10^(snr/20);rng(seed,'twister');noisy.V_cap=daq.V_clean+sigma*randn(size(daq.V_clean));
else
    sigma=0;noisy.V_cap=daq.V_clean;
end
map=map_highspeed_to_space(noisy,cfg.alpha_k,cfg.R_tip,lib.domain,.02);
meta=struct('seed',seed,'signal_rms_V',signal,'noise_std_V',sigma);
end

function [best,meta]=oracle_pair(map,model,cfg,fullFit,pair,highBudget)
q=cfg;q.dualSyncCandidateEOPairs=pair;q.funnelFullRefineCount=1;
if highBudget,q.dualSyncMaxIter=600;else,q.dualSyncMaxIter=240;end
fixed=run_dual_sync_voltage_vp_funnel(map,model,q,fullFit.staticState);
zRef=fixed.fit.theta(:).';fRot=cfg.RPM_high/60;
if ~highBudget
    best=fixed.fit;
    meta=struct('warm_start_rmse_V',fixed.rmse,'oracle_rmse_V',best.rmse,'improvement_V',0);
    return;
end
gLo=max(.05,fullFit.staticState.gHat-cfg.route30GapHalfWidthMm);
gHi=fullFit.staticState.gHat+cfg.route30GapHalfWidthMm;
lb=[gLo,-.5,-1.5,-1.5,-1.5,-1.5];ub=[gHi,.5,1.5,1.5,1.5,1.5];
rng(830001+100*pair(1)+pair(2),'twister');starts=[zRef;...
    repmat(zRef,3,1)+[.08*randn(3,1),.12*randn(3,1),.10*randn(3,4)]];
maxIter=600;
starts=min(max(starts,lb),ub);
% The staged fixed-pair result is itself a valid oracle start.  A random
% restart is accepted only when it improves that already converged basin.
best=fixed.fit;
for i=1:size(starts,1)
    [z,info]=solve_lsq_bounded(@(v)pair_residual(v,pair,map,model),starts(i,:),lb,ub,maxIter,1e-12,1e-12);
    r=pair_residual(z,pair,map,model);rmse=sqrt(mean(r.^2));
    if rmse<best.rmse
        best=unpack_oracle(z,pair,fRot,rmse,info);
    end
end
meta=struct('warm_start_rmse_V',fixed.rmse,'oracle_rmse_V',best.rmse,...
    'improvement_V',fixed.rmse-best.rmse);
end

function r=pair_residual(z,pair,map,model)
x=map.x_v(:);V=map.V_a(:);
if isfield(map,'theta_v')&&numel(map.theta_v)==numel(x),theta=map.theta_v(:);else,theta=2*pi*map.t_v(:);end
u=z(3)*sin(pair(1)*theta)+z(4)*cos(pair(1)*theta)+z(5)*sin(pair(2)*theta)+z(6)*cos(pair(2)*theta);
if isfield(map,'S_v'),sid=map.S_v(:);else,sid=[];end
r=V-eval_gap_template(model,z(1),x-z(2)-u,sid);r(~isfinite(r))=10*max(std(V),1e-3);
end

function fit=unpack_oracle(z,pair,fRot,rmse,info)
fit=struct('rmse',rmse,'theta',z,'eo',pair,'g',z(1),'dx',z(2),...
    'A',[hypot(z(3),z(4)) hypot(z(5),z(6))],...
    'phi',[atan2(z(4),z(3)) atan2(z(6),z(5))],...
    'f',pair*fRot,'solve_info',info);
end

function r=pack_row(c,snr,mode,scope,lowIndex,highIndex,low,high,fit,trueFit,wrongFit,trueMeta,wrongMeta,subspace)
r=struct('scope',scope,'condition_id',c.condition_id,'snr_db',snr,'noise_mode',mode,...
    'low_index',lowIndex,'high_index',highIndex,'low_seed',low.seed,'high_seed',high.seed,...
    'low_noise_std_V',low.noise_std_V,'high_noise_std_V',high.noise_std_V,...
    'g0_est_mm',low.g0_est_mm,'low_selected_window_mm',low.selected_window_mm,...
    'winner_eo1',fit.eo_id(1),'winner_eo2',fit.eo_id(2),'winner_rmse_V',fit.rmse,...
    'true_eo1',c.eo1_true,'true_eo2',c.eo2_true,'true_rmse_V',trueFit.rmse,...
    'wrong_eo1',wrongFit.eo(1),'wrong_eo2',wrongFit.eo(2),'wrong_rmse_V',wrongFit.rmse,...
    'wrong_minus_true_rmse_V',wrongFit.rmse-trueFit.rmse,...
    'true_A1_mm',trueFit.A(1),'true_A2_mm',trueFit.A(2),'true_g_mm',trueFit.g,...
    'wrong_A1_mm',wrongFit.A(1),'wrong_A2_mm',wrongFit.A(2),'wrong_g_mm',wrongFit.g,...
    'winner_is_true',all(fit.eo_id==[c.eo1_true c.eo2_true]),...
    'true_warm_improvement_V',trueMeta.improvement_V,'wrong_warm_improvement_V',wrongMeta.improvement_V,...
    'subspace_rho_max',subspace.rho_max,'subspace_min_angle_deg',subspace.min_angle_deg);
end

function M=subspace_metrics(map,model,trueFit,wrongFit)
theta=map.theta_v(:);x=map.x_v(:);
if isfield(map,'S_v'),sid=map.S_v(:);else,sid=[];end
Xt=local_design(theta,x,sid,model,trueFit);Xw=local_design(theta,x,sid,model,wrongFit);
[Qt,~]=qr(Xt,0);[Qw,~]=qr(Xw,0);s=svd(Qt'*Qw);s=min(max(s,0),1);
M=struct('rho_max',max(s),'min_angle_deg',acosd(max(s)));
end

function X=local_design(theta,x,sid,model,fit)
z=fit.theta;u=z(3)*sin(fit.eo(1)*theta)+z(4)*cos(fit.eo(1)*theta)+...
    z(5)*sin(fit.eo(2)*theta)+z(6)*cos(fit.eo(2)*theta);
fx=eval_gap_derivative(model,z(1),x-z(2)-u,sid);
X=-fx.*[sin(fit.eo(1)*theta),cos(fit.eo(1)*theta),sin(fit.eo(2)*theta),cos(fit.eo(2)*theta)];
X(~isfinite(X))=0;
end

function S=summarize_detail(D)
[G,scope,id,snr,mode]=findgroups(D.scope,D.condition_id,D.snr_db,D.noise_mode);
n=max(G);rows=repmat(summary_row(),n,1);
for k=1:n
    q=D(G==k,:);rows(k).scope=scope(k);rows(k).condition_id=id(k);rows(k).snr_db=snr(k);rows(k).noise_mode=mode(k);rows(k).n=height(q);
    rows(k).winner_true_rate=mean(q.winner_is_true);rows(k).true_beats_wrong_rate=mean(q.wrong_minus_true_rmse_V>0);
    rows(k).median_wrong_minus_true_rmse_V=median(q.wrong_minus_true_rmse_V);
    rows(k).mean_wrong_minus_true_rmse_V=mean(q.wrong_minus_true_rmse_V);
    rows(k).low_template_g0_sd_mm=std(q.g0_est_mm);rows(k).dominant_winner=mode_pair(q.winner_eo1,q.winner_eo2);
end
S=struct2table(rows);
end

function r=summary_row()
r=struct('scope',"",'condition_id',"",'snr_db',NaN,'noise_mode',"",'n',0,...
    'winner_true_rate',NaN,'true_beats_wrong_rate',NaN,'median_wrong_minus_true_rmse_V',NaN,...
    'mean_wrong_minus_true_rmse_V',NaN,'low_template_g0_sd_mm',NaN,'dominant_winner',"");
end

function p=mode_pair(a,b)
u=unique([a b],'rows');[~,~,ix]=unique([a b],'rows');[~,j]=max(accumarray(ix,1));p=sprintf('(%d,%d)',u(j,1),u(j,2));
end

function T=build_subspace_audit(D)
% The oracle fits give the local nonlinear operating point.  Compare only
% X01, where the candidate pair (10,18) is the dominant observed competitor.
q=D(D.condition_id=="X01_corner_near",:);rows=repmat(subspace_row(),height(q),1);
for i=1:height(q)
    rows(i).scope=q.scope(i);rows(i).snr_db=q.snr_db(i);rows(i).noise_mode=q.noise_mode(i);
    rows(i).low_index=q.low_index(i);rows(i).high_index=q.high_index(i);
    rows(i).wrong_minus_true_rmse_V=q.wrong_minus_true_rmse_V(i);
    rows(i).rho_max=q.subspace_rho_max(i);rows(i).min_angle_deg=q.subspace_min_angle_deg(i);
    rows(i).interpretation=ternary(q.wrong_minus_true_rmse_V(i)>0,"true lower RSS","wrong lower RSS");
end
T=struct2table(rows);
end

function r=subspace_row()
r=struct('scope',"",'snr_db',NaN,'noise_mode',"",'low_index',NaN,'high_index',NaN,...
    'wrong_minus_true_rmse_V',NaN,'rho_max',NaN,'min_angle_deg',NaN,'interpretation',"");
end

function write_report(~,~,outDir)
fid=fopen(fullfile(outDir,'README.md'),'w');
fprintf(fid,'# Final EO Competition Mechanism Audit\n\n');
fprintf(fid,'Funnel V2 is unchanged. `detail.csv` separates clean/noisy low and high records, and compares the true X01 pair `(10,12)` with the dominant observed competitor `(10,18)` after high-budget multi-start optimization.\n\n');
fprintf(fid,'Positive `wrong_minus_true_rmse_V` means the true pair has lower oracle residual. The crossed X01 matrix contains 5 low-template seeds by 10 high-record seeds for each DAQ SNR.\n\n');
fprintf(fid,'`subspace.csv` reports the largest canonical correlation and smallest principal angle between the true and `(10,18)` local EO sensitivity subspaces. Values of `rho_max` near one or small angles indicate a geometry-level near-equivalence.\n');
fclose(fid);
end

function C=make_conditions()
id=["C01_baseline";"C02_eo_10_12";"X01_corner_near";"X02_corner_lowamp"];
seed=[1;2;17;18];eo1=[10;10;10;15];eo2=[26;12;12;18];A1=[.25;.25;.25;.15];A2=[.15;.15;.10;.10];
gL=[.5;.5;.5;.7];gH=[.6;.6;.7;.5];p1=[pi/4;pi/4;0;0];p2=[-pi/3;-pi/3;pi/2;pi];
C=table(seed,id,eo1,eo2,A1,A2,gL,gH,p1,p2,'VariableNames',{'condition_seed_id','condition_id','eo1_true','eo2_true','A1_true_mm','A2_true_mm','g_low_mm','g_high_mm','phi1_true_rad','phi2_true_rad'});
end

function cfg=apply_condition(cfg,c,fRot)
cfg.f_true=[c.eo1_true c.eo2_true]*fRot;cfg.A_true=[c.A1_true_mm c.A2_true_mm];cfg.phi_true=[c.phi1_true_rad c.phi2_true_rad];
end

function model=make_truth_model(lib,cfg,gLow)
T=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg.alpha_k));
cal=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));model=build_path_template_model(lib,T,cal);
end

function D=make_condition_daq(geometry,lib,c,cfg,domain)
t=geometry.t(:);omega=cfg.RPM_high*2*pi/60;xHalf=.52*diff(domain);u=zeros(size(t));
for j=1:numel(cfg.A_true),u=u+cfg.A_true(j)*sin(2*pi*cfg.f_true(j)*t+cfg.phi_true(j));end
base=min(eval_gap_template(lib,c.g_low_mm,lib.xGrid));V=base*ones(size(t));Vstatic=base*ones(size(t));
for rev=1:cfg.NumRevs_high
    for sensor=1:numel(cfg.alpha_k)
        te=geometry.T_opr_truth(rev)+cfg.alpha_k(sensor)/omega;idx=find(abs(t-te)<=xHalf/(omega*cfg.R_tip));if isempty(idx),continue;end
        xNom=omega*cfg.R_tip*(t(idx)-te);q=xNom-u(idx);inside=q>=domain(1)&q<=domain(2);
        if any(inside),V(idx(inside))=eval_gap_template(lib,c.g_high_mm,q(inside),sensor);end
        inside0=xNom>=domain(1)&xNom<=domain(2);if any(inside0),Vstatic(idx(inside0))=eval_gap_template(lib,c.g_high_mm,xNom(inside0),sensor);end
    end
end
D=geometry;D.V_clean=V;D.V_static=Vstatic;D.V_cap=V;
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end

function y=ternary(tf,a,b)
if tf,y=a;else,y=b;end
end
