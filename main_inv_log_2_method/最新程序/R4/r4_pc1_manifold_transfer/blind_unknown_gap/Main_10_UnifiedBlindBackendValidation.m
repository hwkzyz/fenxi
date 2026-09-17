function Result = Main_10_UnifiedBlindBackendValidation(opts)
%MAIN_10_UNIFIEDBLINDBACKENDVALIDATION R4 validation with unified backend.
% True FE gaps generate controlled waveforms only.  The estimator receives
% observation/model/config, never the true gaps or vibration parameters.

if nargin < 1, opts = struct(); end
root = fileparts(fileparts(mfilename('fullpath'))); % R4/r4_pc1_manifold_transfer
routeDir = fileparts(root);                  % 最新程序/R4
mainDir = fileparts(routeDir);               % 最新程序 (project root)
addpath(mainDir,'-begin');
addpath(fullfile(mainDir,'local_func'),'-begin');
addpath(fullfile(mainDir,'R2','r2_profiled_recovery'),'-begin');
addpath(fullfile(mainDir,'R3','r3_heldout_recovery'),'-begin');
assert_latest_formal_path(mainDir);
legacyLocal = fullfile(mainDir,'旧程序_20260916','identifiability_mechanism_analysis','local_func');
if exist(legacyLocal,'dir') == 7, addpath(legacyLocal,'-begin'); end
opts = defaults(opts, root);
if ~exist(opts.outputDir,'dir'), mkdir(opts.outputDir); end

ctx = load_inv_log_2_project_context(mainDir);
trust = load_fixed_trust_domain(mainDir);
modelDef = get_inv_log_2_model_def();
gapFE = (0.5:0.2:1.5).';
assert(all(ismembertol([opts.cases.gLow], gapFE, 1e-12)) && ...
    all(ismembertol([opts.cases.gHigh], gapFE, 1e-12)), ...
    'All controlled truth gaps must be existing FE gap nodes.');

D = load(fullfile(root,'output','waveform_families','TiltGap_WaveformSurfaceData.mat'),'Data');
D = D.Data; x = D.x_mm(:); gap = D.gap_mm(:); tilt = D.tilt_deg(:); V = D.response_V;
rows = repmat(result_row(),0,1); lowRows = repmat(low_row(),0,1); audit = struct([]);
for ic = 1:numel(opts.cases)
    c = opts.cases(ic);
    targetIdx = find(abs(tilt-c.tilt) < 1e-12, 1);
    if isempty(targetIdx), error('Controlled tilt %.3g deg is not an FE surface.',c.tilt); end
    targetLib = load_target_library(ctx, trust, modelDef, c.tilt, gapFE, opts.gridN);
    manifold = build_leave_one_out_manifold(V,x,gap,targetIdx,opts.gridN,trust.domain);
    truth = generate_truth(targetLib,c,ctx.cfgAna);
    observation = make_observation(truth, targetLib, ctx.cfgAna);
    config = estimator_config(manifold, opts, ctx.cfgAna);
    assert_no_truth_leak(observation, config);
    [lowCandidates, lowProfile] = profile_low_speed(observation.low, manifold, config);
    for k = 1:numel(lowCandidates)
        lowRows(end+1) = pack_low_row(string(c.id), k, lowCandidates(k)); %#ok<AGROW>
    end
    joint = fit_high_speed_candidates(observation, manifold, lowCandidates, config);
    best = joint(1); secondObjective = best.secondHighObjective;
    if numel(joint) > 1, secondObjective = min(secondObjective,joint(2).jointObjective); end
    margin = (secondObjective-best.jointObjective) / max(best.jointObjective,eps);
    status = classify_result(lowCandidates,best,margin,config);
    rows(end+1) = pack_result(c, best, margin, status, lowProfile, truth); %#ok<AGROW>
    audit(ic).case_id = string(c.id); %#ok<AGROW>
    audit(ic).estimator_input_has_truth = has_truth_field(observation);
    audit(ic).gap_bounds_mm = config.gapBounds;
    audit(ic).low_candidate_count = numel(lowCandidates);
end
T = struct2table(rows); L = struct2table(lowRows);
writetable(T,fullfile(opts.outputDir,'case_results.csv'));
writetable(L,fullfile(opts.outputDir,'low_profile_candidates.csv'));
save(fullfile(opts.outputDir,'blind_unknown_gap_results.mat'), ...
    'T','L','audit','opts','-v7.3');
writetable(struct2table(audit),fullfile(opts.outputDir,'interface_audit.csv'));
Result = struct('cases',T,'low_candidates',L,'audit',audit,'output_dir',string(opts.outputDir));
disp(T);
end

function opts = defaults(opts, root)
if ~isfield(opts,'gridN'), opts.gridN = 401; end
if ~isfield(opts,'gapStep'), opts.gapStep = .05; end
if ~isfield(opts,'zCount'), opts.zCount = 61; end
if ~isfield(opts,'dxGrid'), opts.dxGrid = -.10:.02:.10; end
if ~isfield(opts,'maxLowCandidates'), opts.maxLowCandidates = 3; end
if ~isfield(opts,'maxGridStarts'), opts.maxGridStarts = 2500; end
if ~isfield(opts,'eoGrid'), opts.eoGrid = 1:30; end
if ~isfield(opts,'maxEoPairs'), opts.maxEoPairs = 60; end
if ~isfield(opts,'outputDir'), opts.outputDir = fullfile(root,'blind_unknown_gap','output','gaponly_vp_unknown_eo'); end
if ~isfield(opts,'cases')
    opts.cases = [ ...
        struct('id','increase_09_to_11_eo10','tilt',2.0,'gLow',.9,'gHigh',1.1,'eo',10,'A',.25,'phi',pi/4), ...
        struct('id','decrease_11_to_09_eo17','tilt',2.5,'gLow',1.1,'gHigh',.9,'eo',17,'A',.20,'phi',-pi/3), ...
        struct('id','increase_09_to_11_eo7','tilt',1.5,'gLow',.9,'gHigh',1.1,'eo',7,'A',.22,'phi',pi/6), ...
        struct('id','decrease_11_to_09_eo23','tilt',3.0,'gLow',1.1,'gHigh',.9,'eo',23,'A',.18,'phi',-pi/4)];
end
end

function lib = load_target_library(ctx,trust,modelDef,angle,gapFE,gridN)
name = sprintf('直叶片2mm_不同间隙0.5_0.2_1.5_倾斜角%g.txt',angle);
paths = {fullfile(ctx.rootDir,'data','间隙的影响',name), fullfile(fileparts(mfilename('fullpath')),'..','input_fe_tilt',name)};
for k = 1:numel(paths)
    if exist(paths{k},'file') ~= 2, continue; end
    try
        [g,xCell,yCell] = load_stacked_gap_curves(paths{k},gapFE);
        if numel(g) == numel(gapFE) && all(cellfun(@numel,xCell)>1), break; end
    catch
        xCell = []; yCell = [];
    end
end
if isempty(xCell), error('No readable FE file for %.3g deg.',angle); end
lib = make_response_template_library(gapFE,xCell,yCell,NaN,gridN,modelDef,struct('enable',false));
lib = prepare_gap_template_library(restrict_gap_template_domain(lib,trust.domain));
end

function M = build_leave_one_out_manifold(V,x,gap,targetIdx,gridN,domain)
gref = .8; X = [ones(numel(gap),1), 1./gap, log(gap/gref)];
train = setdiff(1:size(V,3),targetIdx); nX = numel(x); Q = zeros(3*nX,numel(train)); S = zeros(nX*numel(gap),numel(train));
for k = 1:numel(train)
    face = V(:,:,train(k));
    B = (X\face.').'; Q(:,k) = B(:); S(:,k) = face(:);
end
mu = mean(S,2); [U,~,~] = svd(S-mu,'econ'); z = U(:,1)'*(S-mu);
[z,zOrder] = sort(z,'ascend'); M = struct('zTrain',z,'QTrain',Q(:,zOrder),'x',x,'gap',gap,'X',X,'gref',gref,'gridN',gridN,'domain',domain);
dz = max(diff(z)); M.zBounds = [z(1)-.5*dz,z(end)+.5*dz];
end

function lib = library_at_z(M,z)
q = interp1(M.zTrain,M.QTrain.',z,'linear','extrap').'; B = reshape(q,numel(M.x),3); Y = B*M.X.';
xCell = repmat({M.x},numel(M.gap),1); yCell = cell(numel(M.gap),1);
for j = 1:numel(M.gap), yCell{j} = Y(:,j); end
lib = make_response_template_library(M.gap,xCell,yCell,NaN,M.gridN,get_inv_log_2_model_def(),struct('enable',false));
lib = prepare_gap_template_library(restrict_gap_template_domain(lib,M.domain));
end

function truth = generate_truth(targetLib,c,cfg)
% This is the only function in the program that uses controlled truth gaps.
cfg = make_cfg(cfg,c); RPM=cfg.RPM_high; omega=RPM*2*pi/60; T=2*pi/omega; dt=1/cfg.fs;
t=(0:dt:(cfg.NumRevs_high+1)*T)'; T_opr=.08*T+(0:cfg.NumRevs_high)'*T;
% Generate with the same synchronous phase convention used by the blind
% estimator.  The operation-time offset is absorbed into the known phase
% origin; it is not supplied as a vibration-parameter prior.
thetaTruth=omega*(t-T_opr(1));
u=c.A*sin(c.eo*thetaTruth+c.phi);
truthPath=make_low_increment_path_model(targetLib,c.gLow);
baseline=min(truthPath.templateLow(:)); V=baseline*ones(size(t)); xHalf=.52*diff(targetLib.domain);
for rev=1:cfg.NumRevs_high
    for sensor=1:numel(cfg.alpha_k)
        te=T_opr(rev)+cfg.alpha_k(sensor)/omega; idx=find(abs(t-te)<=xHalf/(omega*cfg.R_tip));
        xNom=omega*cfg.R_tip*(t(idx)-te); xPhys=xNom-u(idx); keep=xPhys>=targetLib.domain(1)&xPhys<=targetLib.domain(2);
        V(idx(keep))=eval_path_increment_template(truthPath,c.gHigh-c.gLow,xPhys(keep));
    end
end
truth = struct('gLow',c.gLow,'gHigh',c.gHigh,'deltaG',c.gHigh-c.gLow,'lowWave',truthPath.templateLow(:), ...
    't',t,'V',V,'T_opr',T_opr,'eo',c.eo,'A',c.A,'phi',c.phi,'cfg',cfg);
end

function observation = make_observation(truth,targetLib,cfg)
% Deliberately omit truth gaps and vibration truth from the estimator input.
D = struct('t',truth.t,'V_cap',truth.V,'T_opr_truth',truth.T_opr);
map = map_highspeed_to_space(D,cfg.alpha_k,cfg.R_tip,targetLib.domain,.02);
observation = struct('low',struct('x',targetLib.xGrid(:),'V',truth.lowWave(:)), ...
    'highMap',map);
end

function config = estimator_config(M,opts,cfgAna)
fRot = cfgAna.RPM_high/60;
config = struct('gapBounds',[min(M.gap),max(M.gap)],'gapGrid',min(M.gap):opts.gapStep:max(M.gap), ...
    'zGrid',linspace(M.zBounds(1),M.zBounds(2),opts.zCount),'dxGrid',opts.dxGrid, ...
    'maxLowCandidates',opts.maxLowCandidates,'maxGridStarts',opts.maxGridStarts,'lowDxBound',max(abs(opts.dxGrid)), ...
    'highDxBound',.20,'ampBound',1.5,'ambiguityMargin',.01, ...
    'eoGrid',opts.eoGrid,'maxEoPairs',opts.maxEoPairs, ...
    'rotHz',fRot,'numVarproCandidates',max(opts.maxEoPairs,10), ...
    'f1Grid',opts.eoGrid(:)*fRot,'f2Grid',opts.eoGrid(:)*fRot, ...
    'vpUniqueUnorderedPairs',true,'vpCandidateDiversityHz',1.0, ...
    'vpUseStagedGrid',false,'vpRefineEnable',false, ...
    'gapInitProjectionFreqs',opts.eoGrid(:)*fRot);
end

function [candidates,profile] = profile_low_speed(low,M,config)
J = inf(numel(config.zGrid),numel(config.gapGrid),numel(config.dxGrid));
for iz=1:numel(config.zGrid)
    lib=library_at_z(M,config.zGrid(iz));
    for ig=1:numel(config.gapGrid)
        for id=1:numel(config.dxGrid)
            r=low.V-eval_gap_template(lib,config.gapGrid(ig),low.x-config.dxGrid(id)); J(iz,ig,id)=mean(r.^2);
        end
    end
end
[values,order]=sort(J(:),'ascend'); candidates=repmat(struct('z',NaN,'g',NaN,'dx',NaN,'objective',NaN,'rmse',NaN,'rhoZG',NaN),0,1);
for n=1:min(numel(order),config.maxGridStarts)
    [iz,ig,id]=ind2sub(size(J),order(n)); seed=[config.zGrid(iz),config.gapGrid(ig),config.dxGrid(id)];
    if any(arrayfun(@(q) abs(q.z-seed(1))<.04*range(config.zGrid) && abs(q.g-seed(2))<.04 && abs(q.dx-seed(3))<.02,candidates)), continue; end
    f=@(p) low_residual(p,low,M); [p,~]=solve_lsq_bounded(f,seed,[M.zBounds(1),config.gapBounds(1),-config.lowDxBound], ...
        [M.zBounds(2),config.gapBounds(2),config.lowDxBound],300,1e-10,1e-10);
    r=f(p); candidates(end+1)=struct('z',p(1),'g',p(2),'dx',p(3),'objective',mean(r.^2),'rmse',sqrt(mean(r.^2)),'rhoZG',low_sensitivity_corr(p,low,M)); %#ok<AGROW>
    % Candidate separation must be applied after continuous refinement.  A
    % single smooth valley can otherwise appear several times from nearby
    % coarse-grid starts and create a false identifiability warning.
    if numel(candidates) > 1
        current = candidates(end);
        prior = candidates(1:end-1);
        duplicate = arrayfun(@(q) abs(q.z-current.z)<.01*range(M.zBounds) && ...
            abs(q.g-current.g)<.01 && abs(q.dx-current.dx)<.005, prior);
        if any(duplicate), candidates(end) = []; continue; end
    end
    if numel(candidates)>=config.maxLowCandidates, break; end
end
[~,idx]=sort([candidates.objective]); candidates=candidates(idx); profile=struct('minObjective',values(1),'gridObjective',J);
end

function r = low_residual(p,low,M)
lib=library_at_z(M,p(1)); r=low.V-eval_gap_template(lib,p(2),low.x-p(3)); r(~isfinite(r))=10;
end

function rho = low_sensitivity_corr(p,low,M)
dz=max(.001,.002*range(M.zBounds)); dg=.002;
vz=(low_prediction([p(1)+dz,p(2),p(3)],low,M)-low_prediction([p(1)-dz,p(2),p(3)],low,M))/(2*dz);
vg=(low_prediction([p(1),p(2)+dg,p(3)],low,M)-low_prediction([p(1),p(2)-dg,p(3)],low,M))/(2*dg);
rho=corr(vz,vg);
end

function y = low_prediction(p,low,M), lib=library_at_z(M,p(1)); y=eval_gap_template(lib,p(2),low.x-p(3)); end

function joint = fit_high_speed_candidates(observation,M,lowCandidates,config)
joint=repmat(struct('z',NaN,'gLow',NaN,'dxLow',NaN,'gHigh',NaN,'dxHigh',NaN,'A',NaN,'phi',NaN,'deltaF',NaN,'eo',NaN,'lowObjective',NaN,'highObjective',NaN,'jointObjective',NaN,'highRmse',NaN,'secondHighObjective',NaN,'vpDispRmseMm',NaN,'seedVoltageRmseMv',NaN,'backendStatus',""),0,1);
for k=1:numel(lowCandidates)
    lib=library_at_z(M,lowCandidates(k).z);
    lowTemplate=eval_gap_template(lib,lowCandidates(k).g,lib.xGrid-lowCandidates(k).dx);
    pathModel=make_low_increment_path_model(lib,lowCandidates(k).g,lowTemplate);
    hm=observation.highMap;
    bundle=struct('t',hm.t_v,'x',hm.x_v,'V',hm.V_a,'theta',hm.theta_v, ...
        'sensorId',hm.S_v,'turnId',hm.rev_v,'Wvp',ones(size(hm.V_a)), ...
        'Wfull',ones(size(hm.V_a)));
    response=struct('evalF',@(dg,xq)eval_path_increment_template(pathModel,dg,xq), ...
        'evalFx',@(dg,xq)eval_gap_derivative(pathModel,dg,xq));
    deltaBounds=config.gapBounds-lowCandidates(k).g;
    deltaGrid=config.gapGrid-lowCandidates(k).g;
    bc=struct('gapBounds',deltaBounds,'gapGrid',deltaGrid, ...
        'eoGrid',config.eoGrid,'topK',min(5,config.maxEoPairs), ...
        'rotHz',config.rotHz,'minGradient',1e-6,'minSamples',30, ...
        'dxBound',config.highDxBound,'ampCoeffBound',config.ampBound, ...
        'deltaFBound',2,'maxIter',300,'ambiguityMargin',config.ambiguityMargin);
    fit=run_unified_blind_dynamic_backend(bundle,response,bc);
    if ~isfield(fit,'best'), continue; end
    b=fit.best;
    joint(end+1)=struct('z',lowCandidates(k).z,'gLow',lowCandidates(k).g,'dxLow',lowCandidates(k).dx, ...
        'gHigh',lowCandidates(k).g+b.gap,'dxHigh',b.dx,'A',b.A,'phi',b.phi,'deltaF',b.deltaF,'eo',b.eo,'lowObjective',lowCandidates(k).objective, ...
        'highObjective',b.fullWaveSseV2,'jointObjective',lowCandidates(k).objective+b.fullWaveSseV2,'highRmse',b.fullWaveRmseV, ...
        'secondHighObjective',second_objective(fit),'vpDispRmseMm',b.vpDispRmseMm, ...
        'seedVoltageRmseMv',b.seedVoltageRmseMv,'backendStatus',fit.status); %#ok<AGROW>
end
[~,idx]=sort([joint.jointObjective]); joint=joint(idx);
end

function v=second_objective(fit)
v=NaN;
if isfield(fit,'refined') && numel(fit.refined)>1, v=fit.refined(2).fullWaveSseV2; end
end

function fit = fit_high_speed(map,eoPairs,lib,config)
lb=[config.gapBounds(1),-config.highDxBound,-config.ampBound,-config.ampBound,-config.ampBound,-config.ampBound];
ub=[config.gapBounds(2), config.highDxBound, config.ampBound, config.ampBound, config.ampBound, config.ampBound];
starts=[mean(config.gapBounds),0,0,0,0,0; config.gapBounds(1),0,.1,0,.1,0; config.gapBounds(2),0,-.1,.1,.1,-.1];
best=[]; allObj=inf(size(eoPairs,1),1); for ie=1:size(eoPairs,1)
  eo=eoPairs(ie,:);
  for k=1:size(starts,1)
    [z,si]=solve_lsq_bounded(@(p) high_residual(p,eo,map,lib),starts(k,:),lb,ub,300,1e-10,1e-10); r=high_residual(z,eo,map,lib);
    obj=mean(r.^2); allObj(ie)=min(allObj(ie),obj);
    if isempty(best)||obj<best.objective, best=pack_high_fit(z,eo,r,si); end
  end
end
fit=best; s=sort(allObj(isfinite(allObj))); fit.secondObjective=inf; if numel(s)>1, fit.secondObjective=s(2); end
end

function pairs = screen_eo_pairs(map,lib,config)
% VP screening over admissible static gaps. Frequencies are converted back
% to EO candidates; no true EO or true gap is supplied to this stage.
scores = zeros(0,3);
for ig = 1:numel(config.gapGrid)
    fit = fit_vp_main(map.t_v(:),map.V_a(:),map.x_v(:),lib,config.gapGrid(ig),config);
    for k = 1:numel(fit.candidates)
        eo = round(sort(fit.candidates(k).f(:)/config.rotHz));
        if numel(eo) ~= 2 || any(~ismember(eo,config.eoGrid)) || eo(1)==eo(2), continue; end
        scores(end+1,:) = [fit.candidates(k).linear_sse,eo(:).']; %#ok<AGROW>
    end
end
if isempty(scores)
    pairs = nchoosek(config.eoGrid(:),2);
    return;
end
[~,ord] = sort(scores(:,1),'ascend'); scores = scores(ord,:);
pairs = zeros(0,2);
for k = 1:size(scores,1)
    p = scores(k,2:3);
    if ~ismember(p,pairs,'rows'), pairs(end+1,:) = p; end %#ok<AGROW>
    if size(pairs,1) >= config.maxEoPairs, break; end
end
end

function r = high_residual(p,eo,map,lib)
theta=map.theta_v(:); u=p(3)*sin(eo(1)*theta)+p(4)*cos(eo(1)*theta)+p(5)*sin(eo(2)*theta)+p(6)*cos(eo(2)*theta);
r=map.V_a(:)-eval_gap_template(lib,p(1),map.x_v(:)-p(2)-u); r(~isfinite(r))=10*max(std(map.V_a),1e-3);
end

function fit = pack_high_fit(p,eo,r,si)
fit=struct('g',p(1),'dx',p(2),'eo',eo,'A',[hypot(p(3),p(4)),hypot(p(5),p(6))], ...
    'phi',[atan2(p(4),p(3)),atan2(p(6),p(5))],'rmse',sqrt(mean(r.^2)),'objective',mean(r.^2),'solveInfo',si);
end

function cfg = make_cfg(cfg,c)
cfg.RPM_low=min(cfg.RPM_high,300); cfg.NumRevs_high=8; fRot=cfg.RPM_high/60;
cfg.A_true=c.A; cfg.f_true=c.eo*fRot; cfg.phi_true=c.phi;
end

function status = classify_result(low,best,margin,config)
if abs(low(1).rhoZG)>.98, status="ambiguous_low_speed";
elseif ~isfinite(margin), status="ambiguous_eo";
elseif isfinite(margin)&&margin<config.ambiguityMargin, status="ambiguous_joint";
else, status="identified"; end
end

function r = pack_result(c,best,margin,status,profile,truth)
    r=struct('case_id',string(c.id),'tilt_true_deg',c.tilt,'g_low_true_mm',truth.gLow,'g_low_hat_mm',best.gLow, ...
    'g_low_error_mm',best.gLow-truth.gLow,'g_high_true_mm',truth.gHigh,'g_high_hat_mm',best.gHigh, ...
    'g_high_error_mm',best.gHigh-truth.gHigh,'delta_g_true_mm',truth.deltaG,'delta_g_hat_mm',best.gHigh-best.gLow, ...
    'delta_g_error_mm',best.gHigh-best.gLow-truth.deltaG,'z_hat',best.z,'low_rmse_V',sqrt(best.lowObjective), ...
    'high_rmse_V',best.highRmse,'joint_objective',best.jointObjective,'second_best_margin',margin, ...
    'vp_disp_rmse_mm',best.vpDispRmseMm,'seed_voltage_rmse_mV',best.seedVoltageRmseMv,'delta_f_hat_Hz',best.deltaF, ...
    'eo_true',string(mat2str(c.eo)),'eo_hat',string(mat2str(best.eo)), ...
    'eo_max_abs_error',max(abs(sort(c.eo)-sort(best.eo))), ...
    'A_true',string(mat2str(c.A,5)),'A_hat',string(mat2str(best.A,5)), ...
    'phi_true',string(mat2str(c.phi,5)),'phi_hat',string(mat2str(best.phi,5)), ...
    'low_profile_grid_min',profile.minObjective,'identification_method',"pc1_manifold_vp_fullwave", ...
    'forward_model_mode',"low_increment",'identifiability_status',string(status));
end

function r = result_row(), r=struct('case_id',"",'tilt_true_deg',NaN,'g_low_true_mm',NaN,'g_low_hat_mm',NaN,'g_low_error_mm',NaN,'g_high_true_mm',NaN,'g_high_hat_mm',NaN,'g_high_error_mm',NaN,'delta_g_true_mm',NaN,'delta_g_hat_mm',NaN,'delta_g_error_mm',NaN,'z_hat',NaN,'low_rmse_V',NaN,'high_rmse_V',NaN,'joint_objective',NaN,'second_best_margin',NaN,'vp_disp_rmse_mm',NaN,'seed_voltage_rmse_mV',NaN,'delta_f_hat_Hz',NaN,'eo_true',"",'eo_hat',"",'eo_max_abs_error',NaN,'A_true',"",'A_hat',"",'phi_true',"",'phi_hat',"",'low_profile_grid_min',NaN,'identification_method',"",'forward_model_mode',"",'identifiability_status',""); end
function r = low_row(), r=struct('case_id',"",'rank',NaN,'z_hat',NaN,'g_low_hat_mm',NaN,'dx_low_hat_mm',NaN,'low_objective',NaN,'low_rmse_V',NaN,'rho_zg',NaN); end
function r = pack_low_row(caseId,k,c), r=struct('case_id',caseId,'rank',k,'z_hat',c.z,'g_low_hat_mm',c.g,'dx_low_hat_mm',c.dx,'low_objective',c.objective,'low_rmse_V',c.rmse,'rho_zg',c.rhoZG); end
function assert_no_truth_leak(observation,config), assert(~has_truth_field(observation)&&~has_truth_field(config),'Truth fields leaked into estimator input.'); end
function tf = has_truth_field(s)
tf=false;
if isstruct(s)
    names=string(fieldnames(s));
    if any(contains(names,["truth","gLowTrue","gHighTrue","eoTrue","ATrue","phiTrue"],IgnoreCase=true)), tf=true; return; end
    for k=1:numel(s)
        fn=fieldnames(s(k));
        for j=1:numel(fn)
            v=s(k).(fn{j});
            if isstruct(v) && has_truth_field(v), tf=true; return; end
        end
    end
elseif iscell(s)
    for k=1:numel(s)
        if has_truth_field(s{k}), tf=true; return; end
    end
end
end
