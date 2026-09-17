function R = run_main_vp_joint(highMap, templateLib, cfg, staticState)
%run_main_vp_joint
% Fast static gap init -> VP-style local gap screening ->
% fixed-gap nonlinear refinement -> final projected joint refinement.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);

if nargin < 4 || isempty(staticState) || ~isfield(staticState, 'gHat')
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
end

gCenter = staticState.gHat;
gapHalfWidth = get_cfg_field(cfg, 'mainVpGapHalfWidth', 0.06);
gapGrid = linspace(max(0.05, gCenter - gapHalfWidth), gCenter + gapHalfWidth, ...
    get_cfg_field(cfg, 'mainVpGapN', 61));
if isfield(templateLib, 'fixedPathIncrement') && templateLib.fixedPathIncrement
    gLo=gCenter-gapHalfWidth;gHi=gCenter+gapHalfWidth;
    if isfield(templateLib,'pathCache')&&isfield(templateLib.pathCache,'gGrid')
        % Keep centered finite-difference derivative queries strictly inside
        % the finite cache. Exact cache endpoints can return NaN.
        cacheInset=1e-6;
        gLo=max(gLo,min(templateLib.pathCache.gGrid)+cacheInset);
        gHi=min(gHi,max(templateLib.pathCache.gGrid)-cacheInset);
    end
    gapGrid = linspace(gLo,gHi,get_cfg_field(cfg, 'mainVpGapN', 61));
else
    gapGrid = gapGrid(gapGrid >= min(templateLib.gapTrain) - cfg.rawGapSearchMargin & ...
        gapGrid <= max(templateLib.gapTrain) + cfg.rawGapSearchMargin);
end

timerPrecomp = tic;
precomp = build_vp_precompute(t, V, x, templateLib, gapGrid, cfg);
timePrecomp = toc(timerPrecomp);

bestFit = [];
bestG = gapGrid(1);
bestSSE = inf;
coarseCount = min(get_cfg_field(cfg, 'mainVpCoarseCount', 21), numel(gapGrid));
coarseIdx = round(linspace(1, numel(gapGrid), coarseCount));
coarseIdx = unique(coarseIdx(:));
coarseRows = zeros(numel(coarseIdx), 3);
frequencyPairEvaluations = 0;

timerScreen = tic;
for ic = 1:numel(coarseIdx)
    ig = coarseIdx(ic);
    gTry = gapGrid(ig);
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precomp, ig);
    frequencyPairEvaluations = frequencyPairEvaluations + fitTry.search_pair_count;
    coarseRows(ic, :) = [gTry, fitTry.J1, fitTry.rhoJ];
end

[~, coarseOrder] = sort(coarseRows(:, 2), 'ascend');
bestCoarseJ = coarseRows(coarseOrder(1), 2);
keepCoarse = min(get_cfg_field(cfg, 'mainVpKeepCoarse', 5), numel(coarseOrder));
relTol = get_cfg_field(cfg, 'mainVpRelTol', 0.02);
extraKeep = coarseOrder(coarseRows(coarseOrder, 2) <= bestCoarseJ * (1 + relTol));
keepList = unique([coarseOrder(1:keepCoarse); extraKeep(:)], 'stable');
keepList = keepList(1:min(numel(keepList), get_cfg_field(cfg, 'mainVpMaxKeep', 7)));
fineMask = false(numel(gapGrid), 1);
halfWindow = get_cfg_field(cfg, 'mainVpFineHalfWindow', 5);
for ic = 1:numel(keepList)
    ig0 = coarseIdx(keepList(ic));
    lo = max(1, ig0 - halfWindow);
    hi = min(numel(gapGrid), ig0 + halfWindow);
    fineMask(lo:hi) = true;
end
fineIdx = find(fineMask);

scoreRows = zeros(numel(fineIdx), 3);
fitBank = cell(numel(fineIdx), 1);
for jf = 1:numel(fineIdx)
    ig = fineIdx(jf);
    gTry = gapGrid(ig);
    fitTry = fit_vp_main(t, V, x, templateLib, gTry, cfg, precomp, ig);
    frequencyPairEvaluations = frequencyPairEvaluations + fitTry.search_pair_count;
    scoreRows(jf, :) = [gTry, fitTry.J1, fitTry.rhoJ];
    fitBank{jf} = fitTry;
end
if get_cfg_field(cfg,'mainVpGlobalFreqRefine',false)
    [~,gapOrder]=sort(scoreRows(:,2),'ascend');
    gapOrder=gapOrder(1:min(get_cfg_field(cfg,'mainVpGlobalFreqRefineGapCount',3),numel(gapOrder)));
    for kk=1:numel(gapOrder)
        jf=gapOrder(kk); baseFit=fitBank{jf}; keep=min(4,numel(baseFit.candidates)); f1=[];f2=[];
        for ic=1:keep
            f1=[f1,baseFit.candidates(ic).f(1)+(-20:5:20)]; %#ok<AGROW>
            f2=[f2,baseFit.candidates(ic).f(2)+(-20:5:20)]; %#ok<AGROW>
        end
        cfgLocal=cfg; cfgLocal.vpUseStagedGrid=false; cfgLocal.vpCandidateFrequencyPairsHz=[];
        cfgLocal.f1Grid=unique(f1(f1>=min(cfg.f1Grid)&f1<=max(cfg.f1Grid)));
        cfgLocal.f2Grid=unique(f2(f2>=min(cfg.f2Grid)&f2<=max(cfg.f2Grid)));
        if isempty(cfgLocal.f1Grid)||isempty(cfgLocal.f2Grid),continue;end
        localFit=fit_vp_main(t,V,x,templateLib,scoreRows(jf,1),cfgLocal);
        frequencyPairEvaluations = frequencyPairEvaluations + localFit.search_pair_count;
        if localFit.J1<scoreRows(jf,2),fitBank{jf}=localFit;scoreRows(jf,:)=[scoreRows(jf,1),localFit.J1,localFit.rhoJ];end
    end
end
% Preserve all retained VP candidates for optional full-waveform replay in
% the unified entry. Without this pool, a candidate that ranks poorly on the
% screening subset cannot be recovered later.
candidateBank = repmat(struct('g',NaN,'dx',NaN,'p',[],'f',[],'A',[],'phi',[],...
    'subset_J1',NaN), 0, 1);
for jf = 1:numel(fitBank)
    if isempty(fitBank{jf}) || ~isfield(fitBank{jf}, 'candidates'), continue; end
    localCandidates = fitBank{jf}.candidates;
    for ic = 1:numel(localCandidates)
        c = localCandidates(ic);
        candidateBank(end+1,1) = struct('g',scoreRows(jf,1),'dx',c.dx,...
            'p',c.p,'f',c.f,'A',c.A,'phi',c.phi,...
            'subset_J1',c.linear_sse); %#ok<AGROW>
    end
end

% Conditional completion: dual fits often identify one component reliably
% while the other is replaced by a competing valley. Anchor several strong
% retained frequencies and rescan only the missing component over the full
% one-dimensional grid. This costs O(K*Nf), not another O(Nf^2) search.
if get_cfg_field(cfg, 'mainVpConditionalCompletion', false) && ~isempty(candidateBank)
    [~, baseOrder] = sort([candidateBank.subset_J1], 'ascend');
    anchorSourceCount = min(get_cfg_field(cfg, 'mainVpConditionalAnchorSourceCount', 4), numel(baseOrder));
    anchorPool = [];
    for ia = 1:anchorSourceCount
        anchorPool = [anchorPool, candidateBank(baseOrder(ia)).f]; %#ok<AGROW>
    end
    externalAnchors = get_cfg_field(cfg, 'mainVpExternalAnchorFrequencies', []);
    fallbackAnchors = get_cfg_field(cfg, 'mainVpConditionalAnchorGrid', []);
    anchorPool = unique([externalAnchors(:).', anchorPool, fallbackAnchors(:).'], 'stable');
    anchorCount = min(get_cfg_field(cfg, 'mainVpConditionalAnchorCount', 4), numel(anchorPool));
    anchorPool = anchorPool(1:anchorCount);
    [~, conditionalGapOrder] = sort(scoreRows(:,2), 'ascend');
    conditionalGapCount = min(get_cfg_field(cfg, 'mainVpConditionalGapCount', 3), numel(conditionalGapOrder));
    conditionalGapOrder = conditionalGapOrder(1:conditionalGapCount);
    for igc = 1:numel(conditionalGapOrder)
        jf = conditionalGapOrder(igc);
        for ia = 1:numel(anchorPool)
            cfgConditional = cfg;
            cfgConditional.vpUseStagedGrid = false;
            cfgConditional.vpCandidateFrequencyPairsHz = [];
            cfgConditional.vpUniqueUnorderedPairs = false;
            cfgConditional.f1Grid = cfg.f1Grid;
            cfgConditional.f2Grid = anchorPool(ia);
            cfgConditional.numVarproCandidates = get_cfg_field(cfg, ...
                'mainVpConditionalCandidateCount', 8);
            conditionalFit = fit_vp_main(t, V, x, templateLib, ...
                scoreRows(jf,1), cfgConditional);
            frequencyPairEvaluations = frequencyPairEvaluations + conditionalFit.search_pair_count;
            for ic = 1:numel(conditionalFit.candidates)
                c = conditionalFit.candidates(ic);
                candidateBank(end+1,1) = struct('g',scoreRows(jf,1),...
                    'dx',c.dx,'p',c.p,'f',c.f,'A',c.A,'phi',c.phi,...
                    'subset_J1',c.linear_sse); %#ok<AGROW>
            end
        end
    end
end
timeScreen = toc(timerScreen);

[~, order] = sort(scoreRows(:, 2), 'ascend');
numRefine = min(get_cfg_field(cfg, 'mainVpRefineCount', 5), numel(order));
timerRefineWave = tic;
for ir = 1:numRefine
    jf = order(ir);
    gTry = scoreRows(jf, 1);
    fitTry = fitBank{jf};
    fitTry = refine_vp_main_fit(fitTry, t, V, x, templateLib, gTry, cfg);
    resTry = V - eval_gap_template(templateLib, gTry, x - fitTry.delta_sample);
    sseTry = dot(resTry, resTry);
    if sseTry < bestSSE
        bestSSE = sseTry;
        bestFit = fitTry;
        bestG = gTry;
    end
end
timeRefineWave = toc(timerRefineWave);

if isempty(bestFit)
    jf = order(1);
    gTry = scoreRows(jf, 1);
    fitTry = fitBank{jf};
    bestFit = refine_vp_main_fit(fitTry, t, V, x, templateLib, gTry, cfg);
    bestG = gTry;
end

muList = get_cfg_field(cfg, 'mainUseProjectedJointMuList', [0, 0.15, 0.75]);
bestJoint = [];
bestRmse = inf;
timerJoint = tic;
for imu = 1:numel(muList)
    jointTry = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, ...
        bestG, bestFit.dx, bestFit.p, muList(imu));
    VTry = eval_gap_template(templateLib, jointTry.g, ...
        x - jointTry.dx - fit_u(jointTry.p, t));
    rmseTry = sqrt(mean((V - VTry).^2));
    jointTry.mu = muList(imu);
    jointTry.wave_rmse = rmseTry;
    if rmseTry < bestRmse
        bestRmse = rmseTry;
        bestJoint = jointTry;
    end
end
timeJoint = toc(timerJoint);

VFit = eval_gap_template(templateLib, bestJoint.g, x - bestJoint.dx - fit_u(bestJoint.p, t));

R.method = "main_vp_joint";
R.g_used = bestJoint.g;
R.dx_used = bestJoint.dx;
R.fit = bestJoint;
R.f_id = bestJoint.f;
R.A_id = bestJoint.A;
R.phi_id = bestJoint.phi;
R.VFit = VFit;
R.rmse = bestRmse;
if isfield(templateLib, 'fixedPathIncrement') && templateLib.fixedPathIncrement
    R.eta_g = NaN;
else
    R.eta_g = eta_gap(templateLib, bestJoint.g);
end
R.deltaJ = bestFit.deltaJ;
R.rhoJ = bestFit.rhoJ;
R.J = scoreRows(:, 2);
R.gap_grid = scoreRows(:, 1);
R.candidate_bank = candidateBank;
R.score_table = array2table(scoreRows, 'VariableNames', {'gap_mm','J1','rhoJ'});
R.static_dx0 = staticState.dx0;
R.refined_gap_count = numRefine;
R.coarse_gap_count = numel(coarseIdx);
R.fine_gap_count = numel(fineIdx);
R.time_precomp_s = timePrecomp;
R.time_screen_s = timeScreen;
R.time_refine_wave_s = timeRefineWave;
R.time_joint_s = timeJoint;
R.frequency_pair_evaluations = frequencyPairEvaluations;
R.retained_candidate_count = numel(candidateBank);
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
