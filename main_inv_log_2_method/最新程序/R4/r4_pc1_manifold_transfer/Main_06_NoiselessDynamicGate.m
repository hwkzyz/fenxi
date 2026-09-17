function Result = Main_06_NoiselessDynamicGate(angles)
%MAIN_06_NOISELESSDYNAMICGATE Blind-EO R4 dynamic validation.
% EO is used only by simulate_tilted_truth. The estimator enumerates the
% complete EO pair grid and never receives c.eo or a truth-based initial fit.
% The calibration library remains the zero-tilt library; each FE angle is
% therefore an operating-geometry mismatch relative to calibration.
if nargin<1 || isempty(angles), angles=[0.5 1 1.5 2 2.5 3 3.5]; end
root=fileparts(mfilename('fullpath'));mainDir=fileparts(root);
addpath(mainDir,'-begin');addpath(fullfile(mainDir,'local_func'),'-begin');
addpath(fullfile(mainDir,'identifiability_mechanism_analysis','local_func'));
ctx=load_inv_log_2_project_context(mainDir);trust=load_fixed_trust_domain(mainDir);
gridN=401;modelDef=get_inv_log_2_model_def();
lib=make_fixed_trust_template_library(ctx.gapList,ctx.xCell,ctx.yCell,NaN,gridN,modelDef,trust);
gapFE=(0.5:0.2:1.5).';gLow=.9;gHigh=1.1;delta=gHigh-gLow;
outDir=fullfile(root,'output','noiseless_dynamic_gate');if ~exist(outDir,'dir'),mkdir(outDir);end
angleDir=fullfile(outDir,'per_angle');if ~exist(angleDir,'dir'),mkdir(angleDir);end
rows=repmat(result_row(),0,1);
for ia=1:numel(angles)
    angle=angles(ia);[filePath,sourceKind]=resolve_tilt_file(ctx.rootDir,angle,gapFE);
    [gTilt,xCell,yCell]=load_stacked_gap_curves(filePath,gapFE);
    tilt=make_response_template_library(gTilt,xCell,yCell,NaN,gridN,modelDef,struct('enable',false));
    tilt=restrict_gap_template_domain(tilt,trust.domain);x=lib.xGrid(:);
    yLow=eval_gap_template(tilt,gLow,x);low=struct('xGrid',x,'templateLow',yLow);
    models=cell(3,1);names=["C1_g0_only";"C2_g0_mu";"C3_FE_oracle"];
    masks={[true false false false],[true true false false]};
    for im=1:2
        [p,~]=calibrate_low_speed_path(lib,low,calibration_options(lib,gLow,masks{im}));
        models{im}=build_path_template_model(lib,yLow,p);
    end
    pc=struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1,'b',0,'sensorIds',1:numel(ctx.cfgAna.alpha_k));
    models{3}=build_path_template_model(tilt,yLow,pc);
    cases={struct('id',"null_00",'eo',[10 26],'A',[0 0],'phi',[0 0]),...
        struct('id',"blind_10_26",'eo',[10 26],'A',[.25 .15],'phi',[pi/4 -pi/3]),...
        struct('id',"blind_07_19",'eo',[7 19],'A',[.22 .12],'phi',[pi/6 -pi/5]),...
        struct('id',"blind_13_27",'eo',[13 27],'A',[.20 .10],'phi',[-pi/4 pi/3])};
    for ic=1:numel(cases)
        c=cases{ic};cfg=make_cfg(ctx.cfgAna,c);
        truth=simulate_tilted_truth(tilt,gHigh,c,cfg);
        map=map_highspeed_to_space(truth,cfg.alpha_k,cfg.R_tip,tilt.domain,.02);
        for im=1:3
            fit=blind_pair_fit(map,models{im},cfg);
            if max(c.A)<=eps
                fit.identification_status="not_identifiable_no_vibration";
            end
            rows(end+1)=struct('tilt_deg',angle,'case_id',c.id,'method',names(im),...
                'source_kind',string(sourceKind),'eo_true',string(mat2str(c.eo)),...
                'eo_hat',string(mat2str(fit.eo)),'eo_error',max(abs(sort(c.eo)-sort(fit.eo))),...
                'g_high_true_mm',gHigh,'g_high_hat_mm',fit.g,...
                'gap_error_mm',fit.g-gHigh,'A1_hat_mm',fit.A(1),...
                'A2_hat_mm',fit.A(2),'A1_error_mm',fit.A(1)-c.A(1),...
                'A2_error_mm',fit.A(2)-c.A(2),'rmse_V',fit.rmse,...
                'second_best_margin',fit.second_best_margin,...
                'identification_status',string(fit.identification_status)); %#ok<AGROW>
        end
    end
    angleRows=rows([rows.tilt_deg]==angle);
    AngleTable=struct2table(angleRows);
    writetable(AngleTable,fullfile(angleDir,sprintf('tilt_%03.1fdeg.csv',angle)));
end
T=struct2table(rows);writetable(T,fullfile(outDir,'noiseless_dynamic_gate.csv'));
save(fullfile(outDir,'noiseless_dynamic_gate.mat'),'T','angles','gLow','gHigh','-v7.3');
Result=struct('table',T,'output_dir',string(outDir));disp(T);
end

function cfg=make_cfg(cfg,c)
cfg.RPM_low=min(cfg.RPM_high,300);cfg.NumRevs_high=8;fRot=cfg.RPM_high/60;
cfg.A_true=c.A;cfg.f_true=c.eo*fRot;cfg.phi_true=c.phi;
cfg.route30ForwardModel="low_increment";cfg.route30FrequencyStructureMode="dual_sync_sync";
% Blind search range: independent of the generated EO pair.
cfg.route30DualFrequencyRangeHz=[5 30]*fRot;
cfg.route30GapHalfWidthMm=.20;cfg.dualSyncGapCount=11;cfg.funnelAllPairsFirstJointStep=true;
cfg.funnelSecondJointPairCount=24;cfg.funnelFullRefineCount=6;cfg.dualSyncMaxIter=500;
cfg.dualSyncCandidateEOPairs=[];
cfg.structuredComponentAmplitudeFloorMm=.05;cfg.returnDualSyncDiagnostics=true;
cfg.dualSyncSearchStrategy="funnel_v1";
end

function D=simulate_tilted_truth(tilt,gHigh,c,cfg)
RPM=cfg.RPM_high;nRev=cfg.NumRevs_high;fs=cfg.fs;omega=RPM*2*pi/60;
T=2*pi/omega;dt=1/fs;t=(0:dt:(nRev+1)*T)';oprDelay=.08*T;T_opr=oprDelay+(0:nRev)'*T;
xHalf=.52*diff(tilt.domain);u=zeros(size(t));
for j=1:2,u=u+c.A(j)*sin(2*pi*cfg.f_true(j)*t+c.phi(j));end
baseline=min(tilt.S(:));V=baseline*ones(size(t));Vstatic=V;
for rev=1:nRev
    for sensor=1:numel(cfg.alpha_k)
        te=T_opr(rev)+cfg.alpha_k(sensor)/omega;idx=find(abs(t-te)<=xHalf/(omega*cfg.R_tip));
        if isempty(idx),continue;end
        xNom=omega*cfg.R_tip*(t(idx)-te);xPhys=xNom-u(idx);inside=xPhys>=tilt.domain(1)&xPhys<=tilt.domain(2);
        if any(inside),V(idx(inside))=eval_gap_template(tilt,gHigh,xPhys(inside));end
        inside0=xNom>=tilt.domain(1)&xNom<=tilt.domain(2);
        if any(inside0),Vstatic(idx(inside0))=eval_gap_template(tilt,gHigh,xNom(inside0));end
    end
end
D=struct('t',t,'V_cap',V,'V_clean',V,'V_static',Vstatic,'T_opr',T_opr,...
    'u_truth',u,'noise_std',0,'signal_range',max(V)-min(V),'noise_mode','noiseless');
end

function fit=blind_pair_fit(map,model,cfg)
% Complete blind dual-sync identification. The center gap is the estimated
% low-speed path state; EO pairs are enumerated internally by the funnel.
q=cfg; q.dualSyncCandidateEOPairs=[]; q.funnelAllPairsFirstJointStep=true;
state=struct('gHat',model.pathCal.g0,'g_low_hat',model.pathCal.g0,...
    'dx0',model.pathCal.tau,'calibration_source',"low_speed_path");
R=run_dual_sync_voltage_vp_funnel(map,model,q,state);
fit=R.fit; fit.second_best_margin=second_candidate_margin(R);
if isfinite(fit.second_best_margin) && fit.second_best_margin>=0.01
    fit.identification_status="identified";
else
    fit.identification_status="ambiguous_eo";
end
end

function margin=second_candidate_margin(R)
margin=NaN;
if isfield(R,'diagnostic_refined_rmse') && numel(R.diagnostic_refined_rmse)>=2
    s=sort(R.diagnostic_refined_rmse(:)); margin=(s(2)-s(1))/max(s(1),eps);
elseif isfield(R,'frequency_solution_margin')
    margin=R.frequency_solution_margin;
end
end

function r=pair_residual(z,eo,map,model)
theta=map.theta_v(:);u=z(3)*sin(eo(1)*theta)+z(4)*cos(eo(1)*theta)+...
    z(5)*sin(eo(2)*theta)+z(6)*cos(eo(2)*theta);
r=map.V_a(:)-eval_gap_template(model,z(1),map.x_v(:)-z(2)-u,map.S_v(:));
r(~isfinite(r))=10*max(std(map.V_a),1e-3);
end

function fit=pack_polished(z,eo,map,model,si)
theta=map.theta_v(:);A=[hypot(z(3),z(4)) hypot(z(5),z(6))];
phi=[atan2(z(4),z(3)) atan2(z(6),z(5))];
u=A(1)*sin(eo(1)*theta+phi(1))+A(2)*sin(eo(2)*theta+phi(2));
v=eval_gap_template(model,z(1),map.x_v(:)-z(2)-u,map.S_v(:));
r=map.V_a(:)-v;fit=struct('g',z(1),'dx',z(2),'A',A,'phi',phi,...
    'eo',eo,'theta',z,'rmse',sqrt(mean(r(isfinite(r)).^2)),'solve_info',si);
end

function opts=calibration_options(lib,gLow,mask)
opts=struct('g0Init',gLow,'muInit',0,'tauInit',0,'zetaInit',1,...
    'g0Lower',min(lib.gapTrain),'g0Upper',max(lib.gapTrain),'muLower',-.15,'muUpper',.15,...
    'tauLower',0,'tauUpper',0,'zetaLower',1,'zetaUpper',1,'fitMask',mask,...
    'fitGainOffset',false,'baseline',0,'kappa',1,'xDomain',lib.domain,'maxIter',500);
end

function [filePath,sourceKind]=resolve_tilt_file(rootDir,angle,gapFE)
name=sprintf('直叶片2mm_不同间隙0.5_0.2_1.5_倾斜角%g.txt',angle);
candidates={fullfile(rootDir,'data','间隙的影响',name),...
    fullfile(fileparts(mfilename('fullpath')),'input_fe_tilt',name),...
    fullfile('F:\Program Files\comsol_model\电容数据',name)};
for i=1:numel(candidates)
    if exist(candidates{i},'file')~=2,continue;end
    try
        [g,x,y]=load_stacked_gap_curves(candidates{i},gapFE); %#ok<ASGLU>
        if numel(g)==numel(gapFE)&&all(cellfun(@numel,x)>1)&&all(cellfun(@numel,y)>1)
            filePath=candidates{i};if i==1,sourceKind="workspace";elseif i==2,sourceKind="audit_local_copy";else,sourceKind="legacy_original";end;return;
        end
    catch
    end
end
error('No readable tilted FE file for %.1f deg.',angle);
end

function r=result_row()
r=struct('tilt_deg',NaN,'case_id',"",'method',"",'source_kind',"",...
    'eo_true',"",'eo_hat',"",'eo_error',NaN,'g_high_true_mm',NaN,'g_high_hat_mm',NaN,...
    'gap_error_mm',NaN,'A1_hat_mm',NaN,'A2_hat_mm',NaN,'A1_error_mm',NaN,...
    'A2_error_mm',NaN,'rmse_V',NaN,'second_best_margin',NaN,...
    'identification_status',"");
end
