function Result = Main_04_TransferredSurfaceDynamicValidation(opts)
%MAIN_04_TRANSFERREDSURFACEDYNAMICVALIDATION
% Interface-audited dynamic validation of the PC1/anchor transferred surface.
%
% The program deliberately keeps the three response spaces explicit:
%   (i) absolute static surfaces V(x,g),
%  (ii) the reference-anchored model T_low(x)+kappa*(V(x,g)-V(x,g_low)),
% (iii) the dynamic high-speed observation generated from the target FE face.
% It does not modify Main_02 or Main_06.  The default anchor is g_low=0.9 mm,
% so that the static anchor and the low-speed reference occupy the same gap.
%
% Outputs are written to output/transferred_surface_dynamic_validation/.

if nargin < 1, opts = struct(); end
root = fileparts(mfilename('fullpath'));
mainDir = fileparts(root);
addpath(mainDir,'-begin');
addpath(fullfile(mainDir,'local_func'),'-begin');
addpath(fullfile(mainDir,'identifiability_mechanism_analysis','local_func'));

if ~isfield(opts,'angles'), opts.angles = [0.5 1 1.5 2 2.5 3 3.5]; end
if ~isfield(opts,'anchorGap'), opts.anchorGap = 0.9; end
if ~isfield(opts,'cases'), opts.cases = {'easy_10_26','X01_10_12'}; end

outDir = fullfile(root,'output','transferred_surface_dynamic_validation');
if ~exist(outDir,'dir'), mkdir(outDir); end

ctx = load_inv_log_2_project_context(mainDir);
trust = load_fixed_trust_domain(mainDir);
modelDef = get_inv_log_2_model_def();
gridN = 401;
commonLib = make_fixed_trust_template_library( ...
    ctx.gapList,ctx.xCell,ctx.yCell,NaN,gridN,modelDef,trust);

gapFE = (0.5:0.2:1.5).';
gLow = opts.anchorGap;
if ~ismembertol(gLow,gapFE,1e-12)
    error('anchorGap must be one of the FE gaps [%s].',sprintf(' %.2f',gapFE));
end
anchorIndex = find(abs(gapFE-gLow)<1e-12,1);

% ---- Interface audit before any dynamic result is accepted. ------------
audit = interface_audit(commonLib,ctx,root,gapFE,gLow);
writetable(struct2table(audit),fullfile(outDir,'interface_audit.csv'));
if ~audit.pass
    save(fullfile(outDir,'interface_audit_failed.mat'),'audit','-v7.3');
    error('R4 interface audit failed; dynamic validation was not run.');
end

% Build the static waveform family once.  Main_02 is reproduced exactly,
% including fold-specific centering, PC1 ordering, interpolation and anchor.
dataPath = fullfile(root,'output','waveform_families','TiltGap_WaveformSurfaceData.mat');
D = load(dataPath,'Data'); D = D.Data;
x = D.x_mm(:); gap = D.gap_mm(:); tilt = D.tilt_deg(:); V = D.response_V;
[nX,nGap,nTilt] = size(V);
if numel(gapFE) ~= nGap || max(abs(gap-gapFE)) > 1e-10
    error('Static family gap grid is inconsistent with dynamic FE gap grid.');
end

% Generate one row per target angle and dynamic case for three models.
rows = repmat(result_row(),0,1);
surfaceRows = repmat(surface_row(),0,1);
for it = 1:nTilt
    targetAngle = tilt(it);
    [targetFile,sourceKind] = resolve_tilt_file(ctx.rootDir,targetAngle,gapFE);
    [~,xCellT,yCellT] = load_stacked_gap_curves(targetFile,gapFE);
    targetLib = make_response_template_library(gapFE,xCellT,yCellT,NaN,gridN,modelDef,struct('enable',false));
    targetLib = restrict_gap_template_domain(targetLib,trust.domain);
    targetLib = prepare_gap_template_library(targetLib);
    xDyn = targetLib.xGrid(:);
    yLow = eval_gap_template(targetLib,gLow,xDyn);

    [transferLib,Bhat,latentHat,anchorRmse] = make_transferred_library( ...
        V,x,gap,tilt,it,anchorIndex,gref_value());
    transferLib = align_library_domain(transferLib,trust.domain);
    transferLib = prepare_gap_template_library(transferLib);
    transferLib = resample_library_to_domain(transferLib,targetLib.domain,gridN);
    commonLibDyn = restrict_gap_template_domain(commonLib,targetLib.domain);
    commonLibDyn = prepare_gap_template_library(commonLibDyn);

    % The three models share the same low-speed target trace and g_low path.
    % This is the reference-anchored interface; only the static surface varies.
    pc = struct('g0',gLow,'mu',0,'tau',0,'zeta',1,'kappa',1, ...
        'b',0,'sensorIds',1:numel(ctx.cfgAna.alpha_k));
    mCommon = build_path_template_model(commonLibDyn,yLow,pc);
    mTransfer = build_path_template_model(transferLib,yLow,pc);
    mOracle = build_path_template_model(targetLib,yLow,pc);

    surfaceRows(end+1) = surface_row_from(targetAngle,gLow,latentHat,anchorRmse, ...
        sourceKind,mean(Bhat(:))); %#ok<AGROW>

    for ic = 1:numel(opts.cases)
        c = make_case(opts.cases{ic});
        cfg = make_cfg(ctx.cfgAna,c);
        truth = simulate_tilted_truth(targetLib,1.1,c,cfg);
        map = map_highspeed_to_space(truth,cfg.alpha_k,cfg.R_tip,targetLib.domain,.02);
        models = {mCommon,mTransfer,mOracle};
        names = {'common_reference','anchor_transferred','FE_oracle'};
        for im = 1:3
            % The EO pair is held fixed here.  The purpose of this gate is
            % strictly to propagate a transferred static surface through the
            % R3-compatible reference-anchored forward model, rather than to
            % confound that effect with a competing-E0 search.
            fit = dynamic_pair_fit(map,models{im},cfg,c,1.1);
            rows(end+1) = pack_result(targetAngle,gLow,c.id,names{im},sourceKind, ...
                fit,c,anchorRmse); %#ok<AGROW>
        end
    end
end

T = struct2table(rows);
S = struct2table(surfaceRows);
writetable(T,fullfile(outDir,'transferred_surface_dynamic_results.csv'));
writetable(S,fullfile(outDir,'transferred_surface_records.csv'));
G = groupsummary(T,{'method','case_id'},'mean', ...
    {'gap_error_mm','A1_error_mm','A2_error_mm','fit_rmse_V'});
writetable(G,fullfile(outDir,'transferred_surface_dynamic_summary.csv'));
save(fullfile(outDir,'transferred_surface_dynamic_results.mat'), ...
    'T','S','audit','opts','gapFE','gLow','dataPath','-v7.3');

Result = struct('table',T,'surface',S,'audit',audit,'output_dir',string(outDir));
disp(audit); disp(T);
fprintf('Outputs written to:\n%s\n',outDir);
end

function audit = interface_audit(commonLib,ctx,root,gapFE,gLow)
audit = struct();
audit.pass = true;
audit.common_x_min = min(commonLib.xGrid); audit.common_x_max = max(commonLib.xGrid);
audit.common_gap_min = min(commonLib.gapTrain); audit.common_gap_max = max(commonLib.gapTrain);
audit.anchor_gap_mm = gLow;
audit.n_sensors = numel(ctx.cfgAna.alpha_k);
audit.dynamic_source_files = 0;
for a = [0.5 1 1.5 2 2.5 3 3.5]
    try
        [p,~] = resolve_tilt_file(ctx.rootDir,a,gapFE); %#ok<ASGLU>
        audit.dynamic_source_files = audit.dynamic_source_files + 1;
    catch
        audit.pass = false;
    end
end
audit.has_reference_anchored_evaluator = exist('eval_fixed_path_template','file') == 2;
audit.has_dynamic_mapping = exist('map_highspeed_to_space','file') == 2;
audit.has_gap_evaluator = exist('eval_gap_template','file') == 2;
audit.pass = audit.pass && audit.has_reference_anchored_evaluator && ...
    audit.has_dynamic_mapping && audit.has_gap_evaluator;
end

function [lib,Bhat,latentHat,anchorRmse] = make_transferred_library(V,x,gap,tilt,targetIdx,anchorIdx,gref)
nX = numel(x); nGap = numel(gap); nTilt = numel(tilt);
X = [ones(nGap,1),1./gap,log(gap/gref)];
trainIdx = setdiff(1:nTilt,targetIdx);
C = zeros(3*nX,numel(trainIdx));
Sface = zeros(nX*nGap,numel(trainIdx));
for k = 1:numel(trainIdx)
    B = (X\V(:,:,trainIdx(k)).').';
    C(:,k) = B(:);
    Vk = V(:,:,trainIdx(k));
    Sface(:,k) = Vk(:);
end
mu = mean(Sface,2); [U,~,~] = svd(Sface-mu,'econ');
z = U(:,1)'*(Sface-mu); [zord,ord] = sort(z,'ascend'); Cord = C(:,ord);
dz = max(diff(zord)); tgrid = linspace(zord(1)-0.5*dz,zord(end)+0.5*dz,401);
anchor = V(:,anchorIdx,targetIdx); predAnchor = zeros(nX,numel(tgrid));
for k = 1:numel(tgrid)
    q = interp1(zord,Cord.',tgrid(k),'linear','extrap').';
    B = reshape(q,nX,3); predAnchor(:,k) = B*X(anchorIdx,:).';
end
d = sqrt(mean((predAnchor-anchor).^2,1)); [anchorRmse,ih] = min(d); latentHat=tgrid(ih);
qhat = interp1(zord,Cord.',latentHat,'linear','extrap').'; Bhat=reshape(qhat,nX,3);
Y = Bhat*X.';
xCell = repmat({x(:)},nGap,1); yCell = cell(nGap,1);
for j=1:nGap, yCell{j}=Y(:,j); end
modelDef = get_inv_log_2_model_def();
lib = make_response_template_library(gap,xCell,yCell,NaN,401,modelDef,struct('enable',false));
end

function lib = align_library_domain(lib,domain)
lib = restrict_gap_template_domain(lib,domain);
end

function lib = resample_library_to_domain(lib,domain,gridN)
if abs(lib.domain(1)-domain(1))<1e-10 && abs(lib.domain(2)-domain(2))<1e-10, return; end
% Rebuild the library on its own support; restrict_gap_template_domain has
% already made the support intersection explicit.  A mismatch is fatal.
error('Transferred and target libraries have different spatial domains: [%g,%g] vs [%g,%g].', ...
    lib.domain(1),lib.domain(2),domain(1),domain(2)); %#ok<INUSD>
end

function fit = dynamic_pair_fit(map,model,cfg,c,gHigh)
% Use identical, truth-centred initialization for all three surfaces.  This
% makes FE_oracle a numerical interface check and attributes residual error
% to the static surface, not to funnel-search basin selection.
pair = c.eo;
fRot = cfg.RPM_high/60;
thetaTime = 2*pi*fRot*map.t_v(:);
offset = angle(mean(exp(1i*(map.theta_v(:)-thetaTime))));
phase = c.phi-pair*offset;
z0 = [gHigh 0 c.A(1)*cos(phase(1)) c.A(1)*sin(phase(1)) ...
    c.A(2)*cos(phase(2)) c.A(2)*sin(phase(2))];
lb=[gHigh-.20 -.20 -1.5 -1.5 -1.5 -1.5]; ub=[gHigh+.20 .20 1.5 1.5 1.5 1.5];
res=@(z) pair_residual(z,pair,map,model);
[z,si]=solve_lsq_bounded(res,z0,lb,ub,800,1e-12,1e-12);
fit=pack_fit(z,pair,map,model,si,gHigh);
end

function r = pair_residual(z,eo,map,model)
theta=map.theta_v(:); u=z(3)*sin(eo(1)*theta)+z(4)*cos(eo(1)*theta)+ ...
    z(5)*sin(eo(2)*theta)+z(6)*cos(eo(2)*theta);
r=map.V_a(:)-eval_gap_template(model,z(1),map.x_v(:)-z(2)-u,map.S_v(:));
r(~isfinite(r))=10*max(std(map.V_a),1e-3);
end

function fit=pack_fit(z,eo,map,model,si,gHigh)
theta=map.theta_v(:); A=[hypot(z(3),z(4)),hypot(z(5),z(6))];
phi=[atan2(z(4),z(3)),atan2(z(6),z(5))];
v=eval_gap_template(model,gHigh,map.x_v(:)-z(2)- ...
    A(1)*sin(eo(1)*theta+phi(1))-A(2)*sin(eo(2)*theta+phi(2)),map.S_v(:));
r=map.V_a(:)-v; fit=struct('g',z(1),'A',A,'phi',phi,'rmse',sqrt(mean(r(isfinite(r)).^2)),'solve_info',si);
end

function c=make_case(id)
switch char(id)
    case 'easy_10_26', c=struct('id',string(id),'eo',[10 26],'A',[.25 .15],'phi',[pi/4 -pi/3]);
    case 'X01_10_12', c=struct('id',string(id),'eo',[10 12],'A',[.25 .10],'phi',[0 pi/2]);
    otherwise, error('Unknown dynamic case %s.',id);
end
end

function cfg=make_cfg(cfg,c)
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_high=8; fRot=cfg.RPM_high/60;
cfg.A_true=c.A; cfg.f_true=c.eo*fRot; cfg.phi_true=c.phi;
cfg.route30ForwardModel="low_increment"; cfg.route30FrequencyStructureMode="dual_sync_sync";
cfg.route30DualFrequencyRangeHz=[min(cfg.f_true) max(cfg.f_true)];
cfg.route30GapHalfWidthMm=.20; cfg.dualSyncGapCount=11; cfg.funnelAllPairsFirstJointStep=true;
cfg.funnelSecondJointPairCount=24; cfg.funnelFullRefineCount=1; cfg.dualSyncMaxIter=500;
cfg.structuredComponentAmplitudeFloorMm=.05; cfg.returnDualSyncDiagnostics=true; cfg.dualSyncSearchStrategy="funnel_v1";
end

function D=simulate_tilted_truth(tilt,gHigh,c,cfg)
RPM=cfg.RPM_high; nRev=cfg.NumRevs_high; fs=cfg.fs; omega=RPM*2*pi/60; T=2*pi/omega; dt=1/fs;
t=(0:dt:(nRev+1)*T)'; oprDelay=.08*T; T_opr=oprDelay+(0:nRev)'*T; xHalf=.52*diff(tilt.domain); u=zeros(size(t));
for j=1:2, u=u+c.A(j)*sin(2*pi*cfg.f_true(j)*t+c.phi(j)); end
baseline=min(tilt.S(:)); V=baseline*ones(size(t)); Vstatic=V;
for rev=1:nRev
    for sensor=1:numel(cfg.alpha_k)
        te=T_opr(rev)+cfg.alpha_k(sensor)/omega; idx=find(abs(t-te)<=xHalf/(omega*cfg.R_tip)); if isempty(idx),continue;end
        xNom=omega*cfg.R_tip*(t(idx)-te); xPhys=xNom-u(idx); inside=xPhys>=tilt.domain(1)&xPhys<=tilt.domain(2);
        if any(inside), V(idx(inside))=eval_gap_template(tilt,gHigh,xPhys(inside)); end
        inside0=xNom>=tilt.domain(1)&xNom<=tilt.domain(2); if any(inside0), Vstatic(idx(inside0))=eval_gap_template(tilt,gHigh,xNom(inside0)); end
    end
end
D=struct('t',t,'V_cap',V,'V_clean',V,'V_static',Vstatic,'T_opr_truth',T_opr,'u_truth',u,'noise_std',0,'signal_range',max(V)-min(V),'noise_mode','noiseless');
end

function [filePath,sourceKind]=resolve_tilt_file(rootDir,angle,gapFE)
name=sprintf('直叶片2mm_不同间隙0.5_0.2_1.5_倾斜角%g.txt',angle);
candidates={fullfile(rootDir,'data','间隙的影响',name),fullfile(fileparts(mfilename('fullpath')),'input_fe_tilt',name),fullfile('F:\Program Files\comsol_model\电容数据',name)};
for i=1:numel(candidates)
    if exist(candidates{i},'file')~=2,continue;end
    try
        [g,x,y]=load_stacked_gap_curves(candidates{i},gapFE); %#ok<ASGLU>
        if numel(g)==numel(gapFE)&&all(cellfun(@numel,x)>1)&&all(cellfun(@numel,y)>1)
            filePath=candidates{i}; if i==1,sourceKind="workspace";elseif i==2,sourceKind="audit_local_copy";else,sourceKind="legacy_original";end; return;
        end
    catch
    end
end
error('No readable tilted FE file for %.1f deg.',angle);
end

function r=result_row()
r=struct('tilt_deg',NaN,'anchor_gap_mm',NaN,'case_id',"",'method',"",'source_kind',"",'g_hat_mm',NaN,'gap_error_mm',NaN,'A1_hat_mm',NaN,'A2_hat_mm',NaN,'A1_error_mm',NaN,'A2_error_mm',NaN,'fit_rmse_V',NaN,'anchor_rmse_V',NaN);
end
function r=surface_row()
r=struct('tilt_deg',NaN,'anchor_gap_mm',NaN,'latent_hat',NaN,'anchor_rmse_V',NaN,'source_kind',"",'mean_Bhat',NaN);
end
function r=surface_row_from(t,g,z,e,s,m)
r=surface_row(); r.tilt_deg=t; r.anchor_gap_mm=g; r.latent_hat=z; r.anchor_rmse_V=e; r.source_kind=s; r.mean_Bhat=m;
end
function r=pack_result(t,g,id,name,source,fit,c,e)
r=result_row(); r.tilt_deg=t; r.anchor_gap_mm=g; r.case_id=string(id); r.method=string(name); r.source_kind=source; r.g_hat_mm=fit.g; r.gap_error_mm=fit.g-1.1; r.A1_hat_mm=fit.A(1); r.A2_hat_mm=fit.A(2); r.A1_error_mm=fit.A(1)-c.A(1); r.A2_error_mm=fit.A(2)-c.A(2); r.fit_rmse_V=fit.rmse; r.anchor_rmse_V=e;
end
function g=gref_value(), g=0.8; end
