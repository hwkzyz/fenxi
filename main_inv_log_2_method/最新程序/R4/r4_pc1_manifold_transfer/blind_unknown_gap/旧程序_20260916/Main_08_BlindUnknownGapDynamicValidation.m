function Result = Main_08_BlindUnknownGapDynamicValidation(opts)
%MAIN_08_BLINDUNKNOWNGAPDYNAMICVALIDATION Blind R4 gap-vibration validation.
% True FE gaps generate controlled waveforms only.  The estimator receives
% observation/model/config, never the true gaps or vibration parameters.

if nargin < 1, opts = struct(); end
root = fileparts(fileparts(mfilename('fullpath')));
mainDir = fileparts(root);
addpath(mainDir,'-begin');
addpath(fullfile(mainDir,'local_func'),'-begin');
addpath(fullfile(mainDir,'identifiability_mechanism_analysis','local_func'));
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
    config = estimator_config(manifold, opts);
    assert_no_truth_leak(observation, config);
    [lowCandidates, lowProfile] = profile_low_speed(observation.low, manifold, config);
    for k = 1:numel(lowCandidates)
        lowRows(end+1) = pack_low_row(string(c.id), k, lowCandidates(k)); %#ok<AGROW>
    end
    joint = fit_high_speed_candidates(observation, manifold, lowCandidates, config);
    best = joint(1); secondObjective = NaN;
    if numel(joint) > 1, secondObjective = joint(2).jointObjective; end
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
if ~isfield(opts,'outputDir'), opts.outputDir = fullfile(root,'blind_unknown_gap','output','smoke_fixed_eo_noiseless'); end
if ~isfield(opts,'cases')
    opts.cases = [ ...
        struct('id','increase_09_to_11','tilt',2.0,'gLow',.9,'gHigh',1.1,'eo',[10 26],'A',[.25 .15],'phi',[pi/4 -pi/3]), ...
        struct('id','decrease_11_to_09','tilt',2.5,'gLow',1.1,'gHigh',.9,'eo',[10 26],'A',[.25 .15],'phi',[pi/4 -pi/3])];
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
u=zeros(size(t)); for j=1:2, u=u+c.A(j)*sin(2*pi*cfg.f_true(j)*t+c.phi(j)); end
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

function config = estimator_config(M,opts)
config = struct('gapBounds',[min(M.gap),max(M.gap)],'gapGrid',min(M.gap):opts.gapStep:max(M.gap), ...
    'zGrid',linspace(M.zBounds(1),M.zBounds(2),opts.zCount),'dxGrid',opts.dxGrid, ...
    'maxLowCandidates',opts.maxLowCandidates,'maxGridStarts',opts.maxGridStarts,'lowDxBound',max(abs(opts.dxGrid)), ...
    'highDxBound',.20,'ampBound',1.5,'ambiguityMargin',.01, ...
    'eoGrid',opts.eoGrid,'maxEoPairs',opts.maxEoPairs);
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
joint=repmat(struct('z',NaN,'gLow',NaN,'dxLow',NaN,'gHigh',NaN,'dxHigh',NaN,'A',[],'phi',[],'eo',[],'lowObjective',NaN,'highObjective',NaN,'jointObjective',NaN,'highRmse',NaN),0,1);
for k=1:numel(lowCandidates)
    lib=library_at_z(M,lowCandidates(k).z);
    lowTemplate=eval_gap_template(lib,lowCandidates(k).g,lib.xGrid-lowCandidates(k).dx);
    pathModel=make_low_increment_path_model(lib,lowCandidates(k).g,lowTemplate);
    eoPairs = screen_eo_pairs(observation.highMap,pathModel,config,lowCandidates(k).g);
    fit=fit_high_speed(observation.highMap,eoPairs,pathModel,config,lowCandidates(k).g);
    joint(end+1)=struct('z',lowCandidates(k).z,'gLow',lowCandidates(k).g,'dxLow',lowCandidates(k).dx, ...
        'gHigh',lowCandidates(k).g+fit.g,'dxHigh',fit.dx,'A',fit.A,'phi',fit.phi,'eo',fit.eo,'lowObjective',lowCandidates(k).objective, ...
        'highObjective',fit.objective,'jointObjective',lowCandidates(k).objective+fit.objective,'highRmse',fit.rmse); %#ok<AGROW>
end
[~,idx]=sort([joint.jointObjective]); joint=joint(idx);
end

function fit = fit_high_speed(map,eoPairs,pathModel,config,gLow)
deltaBounds=config.gapBounds-gLow;
lb=[deltaBounds(1),-config.highDxBound,-config.ampBound,-config.ampBound,-config.ampBound,-config.ampBound];
ub=[deltaBounds(2), config.highDxBound, config.ampBound, config.ampBound, config.ampBound, config.ampBound];
starts=[mean(deltaBounds),0,0,0,0,0; deltaBounds(1),0,.1,0,.1,0; deltaBounds(2),0,-.1,.1,.1,-.1];
best=[]; for ie=1:size(eoPairs,1)
  eo=eoPairs(ie,:);
  for k=1:size(starts,1)
    [z,si]=solve_lsq_bounded(@(p) high_residual(p,eo,map,pathModel),starts(k,:),lb,ub,300,1e-10,1e-10); r=high_residual(z,eo,map,pathModel);
    if isempty(best)||mean(r.^2)<best.objective, best=pack_high_fit(z,eo,r,si); end
  end
end
fit=best;
end

function pairs = screen_eo_pairs(map,pathModel,config,gLow)
% Cheap global EO screening: profile the static gap over the full allowed
% gap grid before selecting the EO candidates.  A single mean-gap probe can
% erase the true pair when gap and vibration have comparable amplitudes.
theta=map.theta_v(:); deltaGrid=config.gapGrid-gLow;
scores=[]; grid=config.eoGrid(:).';
for i=1:numel(grid)
 for j=i+1:numel(grid)
   best=Inf;
   for ig=1:numel(deltaGrid)
     yg=map.V_a(:)-eval_path_increment_template(pathModel,deltaGrid(ig),map.x_v(:));
     H=[sin(grid(i)*theta),cos(grid(i)*theta),sin(grid(j)*theta),cos(grid(j)*theta)];
     best=min(best,mean((yg-H*(H\yg)).^2));
   end
   scores(end+1,:)=[best,grid(i),grid(j)]; %#ok<AGROW>
 end
end
[~,ord]=sort(scores(:,1)); n=min(config.maxEoPairs,numel(ord)); pairs=scores(ord(1:n),2:3);
end

function r = high_residual(p,eo,map,pathModel)
theta=map.theta_v(:); u=p(3)*sin(eo(1)*theta)+p(4)*cos(eo(1)*theta)+p(5)*sin(eo(2)*theta)+p(6)*cos(eo(2)*theta);
r=map.V_a(:)-eval_path_increment_template(pathModel,p(1),map.x_v(:)-p(2)-u); r(~isfinite(r))=10*max(std(map.V_a),1e-3);
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
elseif isfinite(margin)&&margin<config.ambiguityMargin, status="ambiguous_joint";
else, status="identified"; end
end

function r = pack_result(c,best,margin,status,profile,truth)
r=struct('case_id',string(c.id),'tilt_true_deg',c.tilt,'g_low_true_mm',truth.gLow,'g_low_hat_mm',best.gLow, ...
    'g_low_error_mm',best.gLow-truth.gLow,'g_high_true_mm',truth.gHigh,'g_high_hat_mm',best.gHigh, ...
    'g_high_error_mm',best.gHigh-truth.gHigh,'delta_g_true_mm',truth.deltaG,'delta_g_hat_mm',best.gHigh-best.gLow, ...
    'delta_g_error_mm',best.gHigh-best.gLow-truth.deltaG,'z_hat',best.z,'low_rmse_V',sqrt(best.lowObjective), ...
    'high_rmse_V',best.highRmse,'joint_objective',best.jointObjective,'second_best_margin',margin, ...
    'eo_true',string(mat2str(c.eo)),'eo_hat',string(mat2str(best.eo)), ...
    'eo_max_abs_error',max(abs(sort(c.eo)-sort(best.eo))), ...
    'low_profile_grid_min',profile.minObjective,'identifiability_status',string(status));
end

function r = result_row(), r=struct('case_id',"",'tilt_true_deg',NaN,'g_low_true_mm',NaN,'g_low_hat_mm',NaN,'g_low_error_mm',NaN,'g_high_true_mm',NaN,'g_high_hat_mm',NaN,'g_high_error_mm',NaN,'delta_g_true_mm',NaN,'delta_g_hat_mm',NaN,'delta_g_error_mm',NaN,'z_hat',NaN,'low_rmse_V',NaN,'high_rmse_V',NaN,'joint_objective',NaN,'second_best_margin',NaN,'eo_true',"",'eo_hat',"",'eo_max_abs_error',NaN,'low_profile_grid_min',NaN,'identifiability_status',""); end
function r = low_row(), r=struct('case_id',"",'rank',NaN,'z_hat',NaN,'g_low_hat_mm',NaN,'dx_low_hat_mm',NaN,'low_objective',NaN,'low_rmse_V',NaN,'rho_zg',NaN); end
function r = pack_low_row(caseId,k,c), r=struct('case_id',caseId,'rank',k,'z_hat',c.z,'g_low_hat_mm',c.g,'dx_low_hat_mm',c.dx,'low_objective',c.objective,'low_rmse_V',c.rmse,'rho_zg',c.rhoZG); end
function assert_no_truth_leak(observation,config), assert(~has_truth_field(observation)&&~has_truth_field(config),'Truth fields leaked into estimator input.'); end
function tf = has_truth_field(s), tf=isstruct(s)&&any(contains(string(fieldnames(s)),["truth","gLowTrue","gHighTrue"],IgnoreCase=true)); end
