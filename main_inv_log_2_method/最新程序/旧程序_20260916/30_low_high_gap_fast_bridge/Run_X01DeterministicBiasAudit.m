function Audit=Run_X01DeterministicBiasAudit()
%RUN_X01DETERMINISTICBIASAUDIT Locate the clean X01 EO-selection bias.
% This is intentionally a deterministic control: no random low/high noise,
% no Funnel budget sweep, and no model-order scan beyond the true and main
% competing EO pairs.
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
outDir=fullfile(root,'output','x01_deterministic_bias_audit');if ~exist(outDir,'dir'),mkdir(outDir);end
ctx=load_inv_log_2_project_context(mainDir);tr=load_fixed_trust_domain(mainDir);
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,...
    ctx.cfgAna.xGridN,get_inv_log_2_model_def(),tr);
cfg=make_cfg(make_inv_log_2_demo_case(ctx,.2).cfgCase);c=x01_case();fRot=cfg.RPM_high/60;
cfg.f_true=[c.eoTrue]*fRot;cfg.A_true=c.A;cfg.phi_true=c.phi;

% Generate the clean DAQ waveform through the same static library used by
% the formal X01 simulation, then map it to the inverse solver coordinates.
truthModel=identity_model(lib,c.gLow,c.gLow,false,cfg);
geometry=simulate_highspeed_from_low_increment(truthModel,.1,cfg,Inf,'fixed_std',1);
daq=make_x01_daq(geometry,lib,c,cfg);map=map_highspeed_to_space(daq,cfg.alpha_k,cfg.R_tip,lib.domain,.02);

% A: the formal low-speed pipeline, including binning, SG and the cached
% linear response surface. B removes only random/template estimation error.
low=simulate_low_speed_template(@(z)eval_gap_template(lib,c.gLow,z),cfg,lib.domain,Inf,901001,'fixed_std');
cal=calibrate_inv_log_2_low_speed(low,lib,cfg);formal=cal.templateModel;
exactCached=identity_model(lib,c.gLow,c.gLow,true,cfg);
exactDirect=identity_model(lib,c.gLow,c.gLow,false,cfg);
models={formal,exactCached,exactDirect};names=["formal_adaptive_template";"exact_low_cached_static";"exact_low_direct_static"];

rows=repmat(row0(),numel(models),1);
for i=1:numel(models)
    m=models{i};rTrue=true_parameter_residual(map,m,c,cfg);
    trueFit=fixed_pair_fit(map,m,cfg,c.eoTrue,c.gHigh);
    wrongFit=fixed_pair_fit(map,m,cfg,[10 18],c.gHigh);
    [rho,angleDeg,proj12,proj18]=bias_geometry(map,m,c,cfg,rTrue);
    rows(i).model=names(i);rows(i).true_parameter_rmse_V=sqrt(mean(rTrue.^2));
    rows(i).true_oracle_rmse_V=trueFit.rmse;rows(i).wrong_oracle_rmse_V=wrongFit.rmse;
    rows(i).wrong_minus_true_oracle_rmse_V=wrongFit.rmse-trueFit.rmse;
    rows(i).rho_max=rho;rows(i).min_angle_deg=angleDeg;
    rows(i).model_error_projection_eo12_V=proj12;rows(i).model_error_projection_eo18_V=proj18;
end
Audit=struct2table(rows);writetable(Audit,fullfile(outDir,'summary.csv'));
save(fullfile(outDir,'audit.mat'),'Audit','map','cal','c','cfg','-v7.3');
write_report(Audit,outDir);disp(Audit);
end

function cfg=make_cfg(cfg)
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_low=20;cfg.NumRevs_high=8;
cfg.route30ForwardModel="low_increment";cfg.route30LowTemplateMethod="adaptive_sg";
cfg.route30LowTemplateOptions=struct('gridSpacingMm',.02,'minBinCount',5,'minCoverage',.90,...
    'numFolds',5,'candidateWindowMm',[.10 .18 .34 .50 .82 1.22 1.62 2.02 2.42]);
cfg.route30FrequencyStructureMode="dual_sync_sync";cfg.route30DualFrequencyRangeHz=[500 1300];
cfg.route30GapHalfWidthMm=.20;cfg.dualSyncGapCount=11;cfg.funnelAllPairsFirstJointStep=true;
cfg.funnelSecondJointPairCount=24;cfg.funnelFullRefineCount=1;cfg.dualSyncMaxIter=600;
cfg.structuredComponentAmplitudeFloorMm=.05;cfg.returnDualSyncDiagnostics=true;cfg.dualSyncSearchStrategy="funnel_v1";
end

function c=x01_case()
c=struct('eoTrue',[10 12],'A',[.25 .10],'phi',[0 pi/2],'gLow',.5,'gHigh',.7,'id',17);
end

function model=identity_model(lib,gLow,g0,useCache,cfg)
T=repmat(eval_gap_template(lib,gLow,lib.xGrid),1,numel(cfg.alpha_k));
pc=struct('g0',g0,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(cfg.alpha_k));
model=build_path_template_model(lib,T,pc);
if ~useCache,model=rmfield(model,'pathCache');end
end

function D=make_x01_daq(geometry,lib,c,cfg)
t=geometry.t(:);omega=cfg.RPM_high*2*pi/60;xHalf=.52*diff(lib.domain);u=zeros(size(t));
for j=1:2,u=u+c.A(j)*sin(2*pi*cfg.f_true(j)*t+c.phi(j));end
base=min(eval_gap_template(lib,c.gLow,lib.xGrid));V=base*ones(size(t));Vstatic=V;
for rev=1:cfg.NumRevs_high
    for sensor=1:numel(cfg.alpha_k)
        te=geometry.T_opr_truth(rev)+cfg.alpha_k(sensor)/omega;idx=find(abs(t-te)<=xHalf/(omega*cfg.R_tip));if isempty(idx),continue;end
        xNom=omega*cfg.R_tip*(t(idx)-te);q=xNom-u(idx);inside=q>=lib.domain(1)&q<=lib.domain(2);
        if any(inside),V(idx(inside))=eval_gap_template(lib,c.gHigh,q(inside),sensor);end
        inside0=xNom>=lib.domain(1)&xNom<=lib.domain(2);
        if any(inside0),Vstatic(idx(inside0))=eval_gap_template(lib,c.gHigh,xNom(inside0),sensor);end
    end
end
D=geometry;D.V_clean=V;D.V_static=Vstatic;D.V_cap=V;
end

function r=true_parameter_residual(map,model,c,cfg)
theta=map.theta_v(:);x=map.x_v(:);sid=map.S_v(:);
% The DAQ generator uses absolute time while the inverse uses OPR angle.
% Convert the generated phases to the inverse angle origin before replay.
thetaTime=2*pi*(cfg.RPM_high/60)*map.t_v(:);offset=angle(mean(exp(1i*(theta-thetaTime))));
phi=c.phi-c.eoTrue*offset;
z=[c.gHigh 0 c.A(1)*cos(phi(1)) c.A(1)*sin(phi(1)) c.A(2)*cos(phi(2)) c.A(2)*sin(phi(2))];
u=z(3)*sin(c.eoTrue(1)*theta)+z(4)*cos(c.eoTrue(1)*theta)+z(5)*sin(c.eoTrue(2)*theta)+z(6)*cos(c.eoTrue(2)*theta);
r=map.V_a(:)-eval_gap_template(model,z(1),x-z(2)-u,sid);r(~isfinite(r))=0;
end

function fit=fixed_pair_fit(map,model,cfg,pair,gHat)
q=cfg;q.dualSyncCandidateEOPairs=pair;q.funnelFullRefineCount=1;
state=struct('gHat',gHat,'g_low_hat',gHat,'dx0',0);
fit=run_dual_sync_voltage_vp_funnel(map,model,q,state).fit;
end

function [rho,angleDeg,p12,p18]=bias_geometry(map,model,c,~,e)
theta=map.theta_v(:);x=map.x_v(:);sid=map.S_v(:);fx=eval_gap_derivative(model,c.gHigh,x,sid);
common=-fx.*[sin(10*theta),cos(10*theta)];
X12=-fx.*[sin(12*theta),cos(12*theta)];X18=-fx.*[sin(18*theta),cos(18*theta)];
valid=all(isfinite([common X12 X18]),2);common=common(valid,:);X12=X12(valid,:);X18=X18(valid,:);e=e(valid);
[Qc,~]=qr(common,0);X12=X12-Qc*(Qc'*X12);X18=X18-Qc*(Qc'*X18);
[Q12,~]=qr(X12,0);[Q18,~]=qr(X18,0);s=svd(Q12'*Q18);rho=max(s);angleDeg=acosd(min(1,max(s)));
p12=norm(Q12*(Q12'*e))/sqrt(numel(e));p18=norm(Q18*(Q18'*e))/sqrt(numel(e));
end

function r=row0()
r=struct('model',"",'true_parameter_rmse_V',NaN,'true_oracle_rmse_V',NaN,'wrong_oracle_rmse_V',NaN,...
    'wrong_minus_true_oracle_rmse_V',NaN,'rho_max',NaN,'min_angle_deg',NaN,...
    'model_error_projection_eo12_V',NaN,'model_error_projection_eo18_V',NaN);
end

function write_report(~,outDir)
fid=fopen(fullfile(outDir,'README.md'),'w');
fprintf(fid,'# X01 Deterministic Bias Audit\n\n');
fprintf(fid,'The three controls isolate: (1) the formal clean low-speed template, (2) an exact low template with the cached linear static surface, and (3) the exact low template with direct PCHIP static evaluation.\n\n');
fprintf(fid,'`true_parameter_rmse_V` is a no-optimization generator-to-inverse replay test after aligning the DAQ time phase with the OPR angle phase. A near-zero result for the direct exact model validates forward consistency. `wrong_minus_true_oracle_rmse_V < 0` means the incorrect `(10,18)` pair wins after fixed-pair full refinement. The reported canonical correlation conditions out the shared EO10 component; an unconditioned comparison is trivially one because both candidates contain EO10.\n');
fclose(fid);
end
