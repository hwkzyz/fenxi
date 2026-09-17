function R = run_low_high_gap_fast_vp_method(highMap, templateLib, cfg, lowState)
%RUN_LOW_HIGH_GAP_FAST_VP_METHOD Formal Route-30 production solver.
%
% The inverse map is not used. A broad unordered dual-frequency VP scan and
% an independent single-frequency VP scan are both refined with the complete
% voltage model. BIC selects model order. The formal path keeps every sample
% in the trusted spatial domain; optional subset screening is diagnostic only.

if nargin < 4 || isempty(lowState) || ~isfield(lowState, 'gHat')
    if ~get_field(cfg, 'allowHighSpeedGapFallback', false)
        error('invlog2:MissingLowSpeedState', ...
            ['The formal main entry requires a low-speed gap state. ' ...
             'Use run_inv_log_2_low_high_main or explicitly set ' ...
             'cfg.allowHighSpeedGapFallback=true for diagnostics.']);
    end
    lowState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
    calibrationSource = "high_speed_fallback";
else
    calibrationSource = get_field(lowState, 'calibration_source', "low_speed_state");
end

% Remove samples that cannot be evaluated by the trusted path model after
% the configured displacement and offset bounds. This is essential for the
% low-speed-template increment model, whose evaluator returns NaN outside
% the trusted spatial support.
originalMap = highMap;
% The absolute template is defined on the full trusted interpolation domain
% and supports extrapolating evaluation used by the validated Route-30 path.
% Query-safe masking is required for the low-speed increment evaluator, whose
% path support is finite and returns NaN outside that support. Applying the
% mask to the absolute template changes the VP objective before candidate
% generation and can discard the samples that made the legacy route robust.
useQuerySafe = isfield(templateLib, 'fixedPathIncrement') && ...
    templateLib.fixedPathIncrement;
if isfield(cfg, 'route30ForceQuerySafe') && cfg.route30ForceQuerySafe
    useQuerySafe = true;
end
if useQuerySafe
    [fitMap, querySafeInfo] = apply_query_safe_window_mask(highMap, templateLib, cfg);
else
    fitMap = highMap;
    querySafeInfo = struct('enabled', false, 'guard_mm', 0, ...
        'domain', templateLib.domain, 'original_count', numel(highMap.t_v), ...
        'safe_count', numel(highMap.t_v), 'removed_count', 0, ...
        'source', 'not_required_for_absolute_template');
end

cfgRoute = apply_route30_defaults(cfg);
snrDb = read_snr_db(cfg);
baseN = cfgRoute.route30BaseSamples;
if isfinite(snrDb) && snrDb <= cfgRoute.route30VeryLowSnrDb
    baseN = max(baseN, cfgRoute.route30LowSnrBaseSamples);
end

fullMap = originalMap;
samplingMode = string(get_field(cfgRoute, 'route30SamplingMode', 'stratified'));
informationSelection = struct('mode',"disabled");
useAllTrusted = get_field(cfgRoute,'route30UseAllTrustedSamples',true);
if useAllTrusted
    % The formal solver never discards a sample that passed the spatial
    % trust-domain check. Information criteria may diagnose conditioning or
    % prioritize candidate work, but they must not redefine the data set.
    workMap = fitMap;
    informationSelection = struct('mode',"all_trusted_samples", ...
        'full_count',numel(fitMap.t_v),'used',numel(fitMap.t_v), ...
        'screening_requested',get_field(cfgRoute,'route30UseInformationSelection',false));
elseif get_field(cfgRoute,'route30UseInformationSelection',false)
    [workMap,informationSelection] = select_information_sufficient_map(...
        fitMap,templateLib,lowState.gHat,cfgRoute);
elseif lower(samplingMode) == "uniform"
    workMap = subset_map_uniform(fitMap, baseN);
else
    workMap = subset_map_stratified(fitMap, baseN);
end
timer = tic;
[dualFit, singleFit, orderInfo] = fit_both_orders(workMap, templateLib, cfgRoute, lowState);
% Re-rank the retained dual VP pool on the complete waveform before local
% refinement. This prevents the screening subset from deciding the final
% basin when its noise realization is unrepresentative.
if orderInfo.model_order == 2
    replayCfg=cfgRoute;
    if isfield(dualFit,'fast_candidate')&&dualFit.fast_candidate.confident
        replayCfg.mainVpConditionalAnchorSourceCount=get_field(cfgRoute,...
            'route30FastConditionalAnchorSourceCount',2);
        replayCfg.mainVpConditionalAnchorCount=get_field(cfgRoute,...
            'route30FastConditionalAnchorCount',6);
        replayCfg.mainVpConditionalGapCount=get_field(cfgRoute,...
            'route30FastConditionalGapCount',2);
        replayCfg.route30IteratedReplayCount=get_field(cfgRoute,...
            'route30FastIteratedReplayCount',30);
        replayCfg.route30FullReplayRefineCount=get_field(cfgRoute,...
            'route30FastFullReplayRefineCount',8);
        replayCfg.route30FullReplayShortMaxIter=get_field(cfgRoute,...
            'route30FastFullReplayShortMaxIter',30);
    end
    timerAugment=tic;
    dualFit = augment_dual_candidates(dualFit, workMap, templateLib, replayCfg, singleFit.f);
    timeAugment=toc(timerAugment);
    timerReplay=tic;
    dualFit = rerank_dual_candidates_full(dualFit, fullMap, workMap, templateLib, replayCfg);
    timeReplay=toc(timerReplay);
    if isfield(dualFit,'pipeline_timing')
        dualFit.pipeline_timing.candidate_augmentation_s=timeAugment;
        dualFit.pipeline_timing.full_replay_refine_s=timeReplay;
    end
end
fallbackUsed = false;
fallbackReason = "none";

if isfield(dualFit,'A_id')&&numel(dualFit.A_id)>=2&&all(isfinite(dualFit.A_id))
    weakRatio = min(dualFit.A_id) / max(max(dualFit.A_id), eps);
else
    weakRatio = NaN;
end
isVeryLowSnr = isfinite(snrDb) && snrDb <= cfgRoute.route30VeryLowSnrDb;
ambiguous = dualFit.rhoJ < cfgRoute.route30RhoFallbackThreshold;
weakDual = orderInfo.model_order == 2 && weakRatio < cfgRoute.route30WeakAmplitudeRatio;
enableWeakFallback = get_field(cfgRoute, 'route30EnableWeakFallback', false) || isVeryLowSnr;
frequencyInfo = estimate_dual_frequency_information(dualFit,workMap,templateLib);
frequencyStdTarget = get_field(cfgRoute,'route30FrequencyStdTargetHz',1.0);
marginTarget = get_field(cfgRoute,'route30ReplayMarginSigmaTarget',2.0);
informationUncertain = (orderInfo.model_order==2||isVeryLowSnr) && ...
    (frequencyInfo.max_frequency_std_hz>frequencyStdTarget || ...
     frequencyInfo.replay_margin_sigma<marginTarget);
needDense = ~useAllTrusted && ((orderInfo.model_order == 2 && ...
    ((enableWeakFallback && weakDual) || ...
     (cfgRoute.route30EnableRhoFallback && ambiguous))) || informationUncertain);

if needDense && numel(fitMap.t_v) > numel(workMap.t_v)
    if informationUncertain
        ratio=max(1,frequencyInfo.max_frequency_std_hz/max(frequencyStdTarget,eps));
        maxGrowth=get_field(cfgRoute,'route30InformationMaxGrowth',4.0);
        denseN=ceil(numel(workMap.t_v)*min(maxGrowth,max(1.8,ratio^2)));
        denseN=min(numel(fitMap.t_v),denseN);
        fallbackReason="frequency_information";
    elseif isVeryLowSnr || ambiguous
        denseN = cfgRoute.route30VeryLowSnrDenseSamples;
        fallbackReason = "low_frequency_margin";
    else
        denseN = cfgRoute.route30WeakComponentSamples;
        fallbackReason = "weak_component";
    end
    if get_field(cfgRoute,'route30UseInformationSelection',false)
        [denseMap,candidateSelectionInfo]=select_candidate_information_map(...
            fitMap,dualFit,templateLib,denseN,cfgRoute);
        informationSelection.refinement=candidateSelectionInfo;
    elseif lower(samplingMode) == "uniform"
        denseMap = subset_map_uniform(fitMap, denseN);
    else
        denseMap = subset_map_stratified(fitMap, denseN);
    end
    if numel(denseMap.t_v) > numel(workMap.t_v)
        [dualFit, singleFit, orderInfo] = fit_both_orders(denseMap, templateLib, cfgRoute, lowState);
        if orderInfo.model_order==2
            dualFit=augment_dual_candidates(dualFit,denseMap,templateLib,cfgRoute,singleFit.f);
            dualFit=rerank_dual_candidates_full(dualFit,fullMap,denseMap,templateLib,cfgRoute);
            frequencyInfo=estimate_dual_frequency_information(dualFit,denseMap,templateLib);
        end
        workMap = denseMap;
        fallbackUsed = true;
    end
end

if orderInfo.model_order == 1
    R = pack_single_result(singleFit, dualFit, orderInfo, lowState, workMap);
else
    R = dualFit;
    R.model_order = 2;
    R.selected_model = "dual";
    R.singleFit = singleFit;
    R.dualFit = dualFit;
    R.bic_single = orderInfo.bic_single;
    R.bic_dual = orderInfo.bic_dual;
    R.delta_bic = orderInfo.delta_bic;
    R.second_frequency_test = orderInfo.second_frequency_test;
end

% Candidate generation may use a stratified subset, but the public result
% must always contain a prediction and RMSE on the complete input waveform.
R.VFit = evaluate_selected_model(R, fullMap, templateLib);
fullResidual = fullMap.V_a(:) - R.VFit(:);
R.rmse = sqrt(mean(fullResidual(isfinite(fullResidual)).^2));
R.fit.wave_rmse = R.rmse;

R.method = "low_high_gap_fast_vp";
R.route = 30;
R.staticState = lowState;
R.calibration_source = calibrationSource;
R.frequency_search_single_hz = [min(cfgRoute.route30SingleGrid), max(cfgRoute.route30SingleGrid)];
R.frequency_search_dual_hz = [min(cfgRoute.f1Grid), max(cfgRoute.f1Grid)];
R.frequency_grid_step_hz = cfgRoute.route30FrequencyStepHz;
if useAllTrusted
    R.base_sample_count = numel(fitMap.t_v);
else
    R.base_sample_count = min(baseN, numel(fitMap.t_v));
end
R.used_sample_count = numel(workMap.t_v);
R.all_trusted_samples_used = useAllTrusted && ...
    numel(workMap.t_v) == numel(fitMap.t_v);
R.fallback_used = fallbackUsed;
R.fallback_reason = fallbackReason;
R.weak_component_ratio = weakRatio;
R.query_safe = querySafeInfo;
R.information_selection = informationSelection;
R.frequency_information = frequencyInfo;
R.dual_search_skipped = isfield(dualFit,'search_skipped')&&dualFit.search_skipped;
R.elapsed_s = toc(timer);
end

function info=estimate_dual_frequency_information(fit,map,lib)
info=struct('frequency_std_hz',[Inf Inf],'max_frequency_std_hz',Inf,...
    'information_condition',Inf,'replay_margin_sigma',0);
if ~isfield(fit,'fit')||~isfield(fit.fit,'p')||numel(fit.fit.p)<6,return;end
t=map.t_v(:);x=map.x_v(:);V=map.V_a(:);p=fit.fit.p;
t0=mean(t);u=fit_u(p,t);z=x-fit.fit.dx-u;
Fx=eval_gap_derivative(lib,fit.g_used,z);h=1e-3;
Fg=(eval_gap_template(lib,fit.g_used+h,z)-eval_gap_template(lib,fit.g_used-h,z))/(2*h);
s1=sin(2*pi*p(3)*t+p(2));c1=cos(2*pi*p(3)*t+p(2));
s2=sin(2*pi*p(6)*t+p(5));c2=cos(2*pi*p(6)*t+p(5));
J=[Fg,-Fx,-Fx.*s1,-Fx.*p(1).*c1,-Fx.*p(1).*c1.*(2*pi*(t-t0)),...
    -Fx.*s2,-Fx.*p(4).*c2,-Fx.*p(4).*c2.*(2*pi*(t-t0))];
V0=eval_gap_template(lib,fit.g_used,z);r=V-V0;
valid=all(isfinite(J),2)&isfinite(r);J=J(valid,:);r=r(valid);
if size(J,1)<=size(J,2),return;end
sigma2=sum(r.^2)/max(size(J,1)-size(J,2),1);
I=J.'*J;C=sigma2*pinv(I);
fstd=sqrt(max(0,[C(5,5),C(8,8)]));
info.frequency_std_hz=fstd;info.max_frequency_std_hz=max(fstd);
info.information_condition=cond(I);
if isfield(fit,'full_replay_scores')&&numel(fit.full_replay_scores)>=2
    s=sort(fit.full_replay_scores(isfinite(fit.full_replay_scores)));
    if numel(s)>=2
        mseStd=max(s(1),eps)*sqrt(2/max(size(J,1)-size(J,2),1));
        info.replay_margin_sigma=(s(2)-s(1))/max(mseStd,eps);
    end
end
end

function VFit = evaluate_selected_model(R, map, lib)
t = map.t_v(:);
x = map.x_v(:);
if R.model_order == 1
    u = R.A_id*sin(2*pi*R.f_id*t+R.phi_id);
else
    u = fit_u(R.fit.p, t);
end
VFit = eval_gap_template(lib, R.g_used, x-R.dx_used-u);
end

function [dualFit, singleFit, info] = fit_both_orders(map, lib, cfg, lowState)
% The verified Route-30 single scan uses one full-band frequency column and
% one fixed nuisance column. Only the leading physical column seeds the
% complete one-frequency voltage refit below.
cfgSingle = cfg;
cfgSingle.f1Grid = cfg.route30SingleGrid;
cfgSingle.f2Grid = cfg.route30SingleNuisanceHz;
cfgSingle.vpUniqueUnorderedPairs = false;
cfgSingle.mainVpGlobalFreqRefine = false;
timerSingle=tic;
singleSeed = fit_single_vp_main(map, lib, cfgSingle, lowState);
singleFit = refit_route30_single_model(map, lib, cfg, singleSeed);
timeSingle=toc(timerSingle);
timerSecond=tic;
secondFrequencyTest = detect_second_frequency_projected(map,singleFit,lib,cfg);
timeSecond=toc(timerSecond);

strongSingle = secondFrequencyTest.delta_bic > ...
    get_field(cfg,'route30SingleEvidenceDeltaBicThreshold',10) && ...
    secondFrequencyTest.amplitude_z < ...
    get_field(cfg,'route30SecondAmplitudeZThreshold',3);
if get_field(cfg,'route30HierarchicalOrderSelection',true)&&strongSingle
    dualFit=struct('search_skipped',true,'skip_reason',"strong_single_evidence",...
        'A_id',[NaN NaN],'deltaJ',NaN,'rhoJ',Inf,...
        'pipeline_timing',struct('single_order_s',timeSingle,...
        'second_frequency_test_s',timeSecond,'candidate_generation_s',0,...
        'dual_gap_profile_s',0,'candidate_augmentation_s',0,...
        'full_replay_refine_s',0));
    n=numel(map.V_a);
    bicSingle=n*log(max(singleFit.rmse^2,eps))+5*log(n);
    info=struct('bic_single',bicSingle,'bic_dual',Inf,'delta_bic',Inf,...
        'model_order',1,'second_frequency_test',secondFrequencyTest);
    return;
end

% The strongest single-frequency component is a cheap, data-driven anchor
% for the conditional dual-frequency completion stage.
cfgDual = cfg;
cfgDual.mainVpExternalAnchorFrequencies = singleFit.f;
cfgDual.mainVpConditionalCompletion = false;
fastInfo=struct('method',"disabled",'confident',false,...
    'fallback_to_full',false,'elapsed_s',0,'candidate_pair_count',0,...
    'full_unordered_pair_count',0,'uses_truth',false);
if get_field(cfg,'route30UseFastDualCandidates',true)
    singleAnchors=singleFit.f;
    if isfield(singleSeed,'candidates')&&~isempty(singleSeed.candidates)
        singleAnchors=[singleAnchors,[singleSeed.candidates.f_id]];
    end
    [fastPairs,fastInfo]=generate_fast_dual_frequency_pairs(...
        map,lib,cfgDual,lowState,singleAnchors);
    if fastInfo.confident
        fastInfo.generated_candidate_pair_count=size(fastPairs,1);
        fastInfo.generated_candidate_pairs_hz=fastPairs;
        fastFreq=unique(fastPairs(:)).';
        cfgDual.f1Grid=fastFreq;cfgDual.f2Grid=fastFreq;
        cfgDual.vpCandidateFrequencyPairsHz=fastPairs;
        cfgDual.mainVpGlobalFreqRefineGapCount=get_field(cfg,...
            'route30FastGlobalRefineGapCount',2);
        cfgDual.mainVpRefineCount=get_field(cfg,'route30FastRefineCount',3);
    end
end
timerDual=tic;
dualFit = run_main_vp_joint(map, lib, cfgDual, lowState);
timeDual=toc(timerDual);
dualFit.fast_candidate=fastInfo;
dualFit.fast_candidate_fallback=fastInfo.fallback_to_full;
dualFit.pipeline_timing=struct('single_order_s',timeSingle,...
    'second_frequency_test_s',timeSecond,...
    'candidate_generation_s',fastInfo.elapsed_s,...
    'dual_gap_profile_s',timeDual,'candidate_augmentation_s',0,...
    'full_replay_refine_s',0);

n = numel(map.V_a);
bicSingle = n * log(max(singleFit.rmse^2, eps)) + 5 * log(n);
bicDual = n * log(max(dualFit.rmse^2, eps)) + 8 * log(n);
info = struct('bic_single', bicSingle, 'bic_dual', bicDual, ...
    'delta_bic', bicDual - bicSingle, ...
    'model_order', 1 + double(bicDual < bicSingle), ...
    'second_frequency_test',secondFrequencyTest);
end

function R = augment_dual_candidates(R0, map, lib, cfg, externalAnchor)
% Add O(K*Nf) conditional candidates to an already computed baseline pool.
R=R0;
if ~get_field(cfg,'mainVpConditionalCompletion',false) || ...
        ~isfield(R0,'candidate_bank') || isempty(R0.candidate_bank)
    return;
end
bank=R0.candidate_bank;
[~,baseOrder]=sort([bank.subset_J1],'ascend');
nSource=min(get_field(cfg,'mainVpConditionalAnchorSourceCount',4),numel(baseOrder));
anchorPool=[];
for i=1:nSource, anchorPool=[anchorPool,bank(baseOrder(i)).f]; end %#ok<AGROW>
anchorPool=unique([externalAnchor(:).',anchorPool,...
    get_field(cfg,'mainVpConditionalAnchorGrid',[])],'stable');
anchorPool=anchorPool(1:min(get_field(cfg,'mainVpConditionalAnchorCount',16),numel(anchorPool)));
[~,gapOrder]=sort(R0.score_table.J1,'ascend');
gapOrder=gapOrder(1:min(get_field(cfg,'mainVpConditionalGapCount',3),numel(gapOrder)));
t=map.t_v(:); V=map.V_a(:); x=map.x_v(:);
% Conditional scans reuse the same full-band basis for every anchor and the
% same forward template for every anchor at a given gap. Cache both so all
% trusted samples still contribute without repeating template/trigonometric
% evaluation K times.
SFull=sin(2*pi*t*cfg.f1Grid(:).');
CFull=cos(2*pi*t*cfg.f1Grid(:).');
SAnchor=sin(2*pi*t*anchorPool);
CAnchor=cos(2*pi*t*anchorPool);
for ig=1:numel(gapOrder)
    g=R0.score_table.gap_mm(gapOrder(ig));
    F0=eval_gap_template(lib,g,x);
    Fx=eval_gap_derivative(lib,g,x);
    precomp=struct('F0Mat',F0,'FxMat',Fx,'YMat',V-F0, ...
        'S1',SFull,'C1',CFull,'S2',[],'C2',[]);
    for ia=1:numel(anchorPool)
        c=cfg; c.vpUseStagedGrid=false; c.vpUniqueUnorderedPairs=false;
        c.f1Grid=cfg.f1Grid; c.f2Grid=anchorPool(ia);
        c.numVarproCandidates=get_field(cfg,'mainVpConditionalCandidateCount',8);
        precomp.S2=SAnchor(:,ia);precomp.C2=CAnchor(:,ia);
        f=fit_vp_main(t,V,x,lib,g,c,precomp,1);
        for k=1:numel(f.candidates)
            q=f.candidates(k);
            bank(end+1,1)=struct('g',g,'dx',q.dx,'p',q.p,'f',q.f,...
                'A',q.A,'phi',q.phi,'subset_J1',q.linear_sse); %#ok<AGROW>
        end
    end
end
R.candidate_bank=bank;
end

function R = rerank_dual_candidates_full(R0, fullMap, workMap, lib, cfg)
%RERANK_DUAL_CANDIDATES_FULL Replay retained VP seeds on all samples, then
% refine only a small number of the best full-waveform seeds.
R = R0;
if ~isfield(R0, 'candidate_bank') || isempty(R0.candidate_bank)
    return;
end
tFull = fullMap.t_v(:); xFull = fullMap.x_v(:); VFull = fullMap.V_a(:);
bank = R0.candidate_bank;
keys = [[bank.g].', vertcat(bank.f)];
[~, uniqueIdx] = unique(round(keys*1e9)/1e9, 'rows', 'stable');
bank = bank(uniqueIdx);
score = inf(numel(bank), 1);
[gapValues,~,gapId] = unique([bank.g].');
F0Cache = cell(numel(gapValues),1);
FxCache = cell(numel(gapValues),1);
for ig = 1:numel(gapValues)
    F0Cache{ig} = eval_gap_template(lib,gapValues(ig),xFull);
    FxCache{ig} = eval_gap_derivative(lib,gapValues(ig),xFull);
end
for i = 1:numel(bank)
    [score(i), bank(i)] = full_linear_replay(bank(i),tFull,VFull,...
        F0Cache{gapId(i)},FxCache{gapId(i)});
end
[~, ord] = sort(score, 'ascend');
% A first-order VP replay can rank large-amplitude candidates poorly. Apply
% a small number of exact-forward Gauss-Newton replay steps to a wider pool
% before spending the expensive joint nonlinear budget.
iterKeep=min(get_field(cfg,'route30IteratedReplayCount',30),numel(ord));
if get_field(cfg,'route30UseIteratedReplay',true)
    for ii=1:iterKeep
        [score(ord(ii)),bank(ord(ii))]=full_replay_iterated(bank(ord(ii)),...
            tFull,xFull,VFull,lib,get_field(cfg,'route30IteratedReplayMaxIter',2));
    end
    [score,ord]=sort(score,'ascend');
end
keep = min(get_field(cfg, 'route30FullReplayRefineCount', 10), numel(ord));
if keep < 1, return; end
t = workMap.t_v(:); x = workMap.x_v(:); V = workMap.V_a(:);
[~,~,gidx] = unique([workMap.rev_v(:), workMap.S_v(:)], 'rows');
ng = max(gidx);
bestRmse = inf;
best = [];
cfgShort = cfg;
cfgShort.projectedJointMaxIter = get_field(cfg,'route30FullReplayShortMaxIter',60);
for i = 1:keep
    c = bank(ord(i));
    try
        joint = refine_projected_joint(t, V, x, lib, cfgShort, gidx, ng, ...
            c.g, c.dx, c.p, 0);
        VFit = eval_gap_template(lib, joint.g, ...
            xFull - joint.dx - fit_u(joint.p, tFull));
        rr = VFull - VFit;
        good = isfinite(rr);
        rmse = sqrt(mean(rr(good).^2));
        if rmse < bestRmse
            bestRmse = rmse;
            best = joint;
        end
    catch
        % Keep the original fit if a replay seed is outside the trust region.
    end
end
if isempty(best), return; end
% Spend the full nonlinear budget only once, on the best short-fit basin.
try
    finalJoint = refine_projected_joint(t,V,x,lib,cfg,gidx,ng,...
        best.g,best.dx,best.p,0);
    finalV = eval_gap_template(lib,finalJoint.g,...
        xFull-finalJoint.dx-fit_u(finalJoint.p,tFull));
    finalResidual=VFull-finalV; good=isfinite(finalResidual);
    finalRmse=sqrt(mean(finalResidual(good).^2));
    if finalRmse <= bestRmse
        best=finalJoint; bestRmse=finalRmse;
    end
catch
    % The best short fit remains valid if the final polish cannot start.
end
R.fit = best;
R.g_used = best.g;
R.dx_used = best.dx;
R.f_id = best.f;
R.A_id = best.A;
R.phi_id = best.phi;
R.VFit = eval_gap_template(lib, best.g, xFull-best.dx-fit_u(best.p,tFull));
R.rmse = bestRmse;
R.full_replay_score = score(ord(1));
R.full_replay_candidate_count = numel(bank);
R.full_replay_scores = score;
R.full_replay_order = ord;
R.full_replay_frequencies = vertcat(bank.f);
end

function [sse, c] = full_linear_replay(c, t, V, F0, Fx)
% Refit dx and sine/cosine coefficients on the complete waveform while
% holding the candidate gap and frequencies fixed.
f=sort(c.f(:).');
S1=sin(2*pi*f(1)*t); C1=cos(2*pi*f(1)*t);
S2=sin(2*pi*f(2)*t); C2=cos(2*pi*f(2)*t);
valid=isfinite(F0)&isfinite(Fx)&isfinite(V);
q=-Fx(valid); y=V(valid)-F0(valid); W=abs(Fx(valid));
X=[q,q.*S1(valid),q.*C1(valid),q.*S2(valid),q.*C2(valid)];
beta=(X.*W)\(y.*W);
r=y-X*beta; sse=mean(r.^2);
if ~isfinite(sse), return; end
A=[hypot(beta(2),beta(3)),hypot(beta(4),beta(5))];
phi=[atan2(beta(3),beta(2)),atan2(beta(5),beta(4))];
c.dx=beta(1); c.f=f; c.A=A; c.phi=phi;
c.p=[A(1),phi(1),f(1),A(2),phi(2),f(2)];
end

function [sse,c]=full_replay_iterated(c,t,x,V,lib,maxIter)
% Exact-forward fixed-frequency Gauss-Newton replay for candidate ranking.
f=sort(c.f(:).');dx=c.dx;coef=[c.A(1)*cos(c.phi(1));c.A(1)*sin(c.phi(1));...
    c.A(2)*cos(c.phi(2));c.A(2)*sin(c.phi(2))];
for it=1:maxIter
    u=coef(1)*sin(2*pi*f(1)*t)+coef(2)*cos(2*pi*f(1)*t)+...
       coef(3)*sin(2*pi*f(2)*t)+coef(4)*cos(2*pi*f(2)*t);
    z=x-dx-u;F0=eval_gap_template(lib,c.g,z);Fx=eval_gap_derivative(lib,c.g,z);
    valid=isfinite(V)&isfinite(F0)&isfinite(Fx);q=-Fx(valid);y=V(valid)-F0(valid);
    S1=sin(2*pi*f(1)*t(valid));C1=cos(2*pi*f(1)*t(valid));
    S2=sin(2*pi*f(2)*t(valid));C2=cos(2*pi*f(2)*t(valid));
    X=[q,q.*S1,q.*C1,q.*S2,q.*C2];delta=X\y;
    dx=dx+delta(1);coef=coef+delta(2:5);
end
u=coef(1)*sin(2*pi*f(1)*t)+coef(2)*cos(2*pi*f(1)*t)+...
   coef(3)*sin(2*pi*f(2)*t)+coef(4)*cos(2*pi*f(2)*t);
VFit=eval_gap_template(lib,c.g,x-dx-u);
r=V-VFit;good=isfinite(r);sse=mean(r(good).^2);
A=[hypot(coef(1),coef(2)),hypot(coef(3),coef(4))];phi=[atan2(coef(2),coef(1)),atan2(coef(4),coef(3))];
c.dx=dx;c.A=A;c.phi=phi;c.f=f;c.p=[A(1),phi(1),f(1),A(2),phi(2),f(2)];
end

function fit = refit_route30_single_model(map, lib, cfg, seed)
t = map.t_v(:);
x = map.x_v(:);
V = map.V_a(:);
[~, component] = max(seed.A_id);
f0 = seed.f_id(component);
a0 = seed.A_id(component) * cos(seed.phi_id(component));
b0 = seed.A_id(component) * sin(seed.phi_id(component));
z0 = [seed.g_used, seed.dx_used, a0, b0, f0];

fMin = min(cfg.route30SingleGrid);
fMax = max(cfg.route30SingleGrid);
freqHalf = cfg.route30SingleRefineHalfWidthHz;
lb = [max(0.05, seed.g_used-cfg.route30SingleGapHalfWidthMm), ...
      -cfg.route30DxHalfWidthMm, -cfg.route30AmplitudeUpperMm, ...
      -cfg.route30AmplitudeUpperMm, max(fMin, f0-freqHalf)];
ub = [seed.g_used+cfg.route30SingleGapHalfWidthMm, ...
      cfg.route30DxHalfWidthMm, cfg.route30AmplitudeUpperMm, ...
      cfg.route30AmplitudeUpperMm, min(fMax, f0+freqHalf)];
if ~(isfield(lib, 'fixedPathIncrement') && lib.fixedPathIncrement)
    lb(1) = max(lb(1), min(lib.gapTrain)-cfg.rawGapSearchMargin);
    ub(1) = min(ub(1), max(lib.gapTrain)+cfg.rawGapSearchMargin);
end
z0 = min(max(z0, lb), ub);

penalty = 10 * max(std(V), 1e-3);
residual = @(z) single_voltage_residual(z, t, V, x, lib, penalty);
[z, solveInfo] = solve_lsq_bounded(residual, z0, lb, ub, ...
    cfg.route30SingleMaxIter, 1e-12, 1e-12);
A = hypot(z(3), z(4));
phi = atan2(z(4), z(3));
VFit = eval_gap_template(lib, z(1), ...
    x-z(2)-A*sin(2*pi*z(5)*t+phi));
r = V - VFit;
valid = isfinite(r);
fit = struct('g', z(1), 'dx', z(2), 'A', A, 'phi', phi, 'f', z(5), ...
    'p', [A, phi, z(5)], 'theta', z, 'VFit', VFit, ...
    'rmse', sqrt(mean(r(valid).^2)), 'solve_info', solveInfo, ...
    'seed_frequency_hz', f0);
end

function r = single_voltage_residual(z, t, V, x, lib, penalty)
A = hypot(z(3), z(4));
phi = atan2(z(4), z(3));
y = eval_gap_template(lib, z(1), x-z(2)-A*sin(2*pi*z(5)*t+phi));
r = V-y;
r(~isfinite(r)) = penalty;
end

function R = pack_single_result(single, dual, info, lowState, map)
R = struct();
R.model_order = 1;
R.selected_model = "single";
R.g_used = single.g;
R.dx_used = single.dx;
R.fit = single;
R.f_id = single.f;
R.A_id = single.A;
R.phi_id = single.phi;
R.VFit = single.VFit;
R.rmse = single.rmse;
R.eta_g = NaN;
R.deltaJ = dual.deltaJ;
R.rhoJ = dual.rhoJ;
R.singleFit = single;
R.dualFit = dual;
R.bic_single = info.bic_single;
R.bic_dual = info.bic_dual;
R.delta_bic = info.delta_bic;
R.second_frequency_test = info.second_frequency_test;
R.static_dx0 = lowState.dx0;
R.sample_time = map.t_v(:);
end

function cfg = apply_route30_defaults(cfg)
freqRange = get_field(cfg, 'route30DualFrequencyRangeHz', [300, 1500]);
singleRange = get_field(cfg, 'route30SingleFrequencyRangeHz', [5, 1500]);
step = get_field(cfg, 'route30FrequencyStepHz', 5);
cfg.route30FrequencyStepHz = step;
cfg.f1Grid = freqRange(1):step:freqRange(2);
cfg.f2Grid = cfg.f1Grid;
cfg.route30SingleGrid = singleRange(1):step:singleRange(2);
cfg.route30SingleNuisanceHz = singleRange(2);
cfg.vpUniqueUnorderedPairs = true;
cfg.vpUseStagedGrid = false;
cfg.numVarproCandidates = get_field(cfg, 'route30CandidateCount', 12);
% The low-speed gap is only a center estimate. A +/-0.70 mm profile covers
% the validated unequal low/high library pairs without assuming equality.
cfg.mainVpGapHalfWidth = get_field(cfg, 'route30GapHalfWidthMm', 0.70);
cfg.mainVpGapN = get_field(cfg, 'route30GapCount', 61);
cfg.mainVpCoarseCount = get_field(cfg, 'route30GapCoarseCount', 31);
cfg.mainVpMaxKeep = get_field(cfg, 'route30GapMaxKeep', 9);
cfg.mainVpGlobalFreqRefine = true;
cfg.mainVpGlobalFreqRefineGapCount = get_field(cfg, 'route30GlobalRefineGapCount', 3);
cfg.mainVpRefineCount = get_field(cfg, 'route30RefineCount', 5);
cfg.mainVpConditionalCompletion = get_field(cfg, 'route30ConditionalCompletion', true);
cfg.mainVpConditionalAnchorSourceCount = get_field(cfg, 'route30ConditionalAnchorSourceCount', 4);
cfg.mainVpConditionalAnchorCount = get_field(cfg, 'route30ConditionalAnchorCount', 16);
cfg.mainVpConditionalGapCount = get_field(cfg, 'route30ConditionalGapCount', 3);
cfg.mainVpConditionalCandidateCount = get_field(cfg, 'route30ConditionalCandidateCount', 8);
cfg.mainVpConditionalAnchorGrid = get_field(cfg, 'route30ConditionalAnchorGridHz', 300:100:1500);
cfg.route30FullReplayRefineCount = get_field(cfg, 'route30FullReplayRefineCount', 10);
cfg.route30FullReplayShortMaxIter = get_field(cfg, 'route30FullReplayShortMaxIter', 60);
cfg.vpCandidateDiversityHz = get_field(cfg, 'route30CandidateDiversityHz', 8);
cfg.mainUseProjectedJointMuList = 0;
cfg.route30BaseSamples = get_field(cfg, 'route30BaseSamples', 900);
cfg.route30SamplingMode = get_field(cfg, 'route30SamplingMode', 'stratified');
cfg.route30UseInformationSelection = get_field(cfg,'route30UseInformationSelection',false);
cfg.route30UseAllTrustedSamples = get_field(cfg,'route30UseAllTrustedSamples',true);
cfg.route30LowSnrBaseSamples = get_field(cfg, 'route30LowSnrBaseSamples', 1500);
cfg.route30VeryLowSnrDenseSamples = get_field(cfg, 'route30VeryLowSnrDenseSamples', 2500);
cfg.route30WeakComponentSamples = get_field(cfg, 'route30WeakComponentSamples', 1500);
cfg.route30VeryLowSnrDb = get_field(cfg, 'route30VeryLowSnrDb', 5);
cfg.route30RhoFallbackThreshold = get_field(cfg, 'route30RhoFallbackThreshold', 0.003);
cfg.route30WeakAmplitudeRatio = get_field(cfg, 'route30WeakAmplitudeRatio', 0.30);
cfg.route30EnableRhoFallback = get_field(cfg, 'route30EnableRhoFallback', true);
cfg.route30EnableWeakFallback = get_field(cfg, 'route30EnableWeakFallback', false);
cfg.route30FrequencyStdTargetHz = get_field(cfg,'route30FrequencyStdTargetHz',1.0);
cfg.route30ReplayMarginSigmaTarget = get_field(cfg,'route30ReplayMarginSigmaTarget',2.0);
cfg.route30InformationMaxGrowth = get_field(cfg,'route30InformationMaxGrowth',4.0);
cfg.route30HierarchicalOrderSelection = get_field(cfg,'route30HierarchicalOrderSelection',true);
cfg.route30SingleEvidenceDeltaBicThreshold = get_field(cfg,...
    'route30SingleEvidenceDeltaBicThreshold',10);
cfg.route30SecondDeltaBicThreshold = get_field(cfg,'route30SecondDeltaBicThreshold',-10);
cfg.route30SecondAmplitudeZThreshold = get_field(cfg,'route30SecondAmplitudeZThreshold',3);
cfg.route30SecondMinSeparationHz = get_field(cfg,'route30SecondMinSeparationHz',20);
cfg.route30SingleRefineHalfWidthHz = get_field(cfg, 'route30SingleRefineHalfWidthHz', 20);
cfg.route30SingleGapHalfWidthMm = get_field(cfg, 'route30SingleGapHalfWidthMm', 0.12);
cfg.route30DxHalfWidthMm = get_field(cfg, 'route30DxHalfWidthMm', 0.5);
cfg.route30AmplitudeUpperMm = get_field(cfg, 'route30AmplitudeUpperMm', 1.5);
cfg.route30SingleMaxIter = get_field(cfg, 'route30SingleMaxIter', 160);
cfg.route30UseFastDualCandidates = get_field(cfg,'route30UseFastDualCandidates',true);
cfg.route30FastCandidateGapAnchors = get_field(cfg,'route30FastCandidateGapAnchors',5);
cfg.route30FastCandidatePairsPerGap = get_field(cfg,'route30FastCandidatePairsPerGap',16);
cfg.route30FastCandidateMaxPairs = get_field(cfg,'route30FastCandidateMaxPairs',640);
end

function map = subset_map_stratified(map, targetN)
n = numel(map.t_v);
if n <= targetN
    map.sampling_info = struct('mode', "complete", 'requested', targetN, ...
        'used', n, 'group_count', count_groups(map));
    return;
end
if ~isfield(map,'rev_v') || ~isfield(map,'S_v')
    idx = unique(round(linspace(1,n,targetN)));
    mode = "uniform_fallback";
    groupCount = NaN;
else
    keys = [map.rev_v(:), map.S_v(:)];
    [groupKeys,~,groupId] = unique(keys,'rows','stable');
    groupCount = size(groupKeys,1);
    minPerGroup = min(4, floor(targetN/max(groupCount,1)));
    minPerGroup = max(minPerGroup,1);
    counts = accumarray(groupId,1,[groupCount,1]);
    raw = targetN * counts / n;
    take = max(minPerGroup, floor(raw));
    take = min(take, counts);
    while sum(take) < min(targetN,n)
        room = counts-take;
        gain = raw-floor(raw);
        gain(room<=0) = -Inf;
        [~,k] = max(gain);
        if ~isfinite(gain(k)), [~,k] = max(room); end
        if room(k)<=0, break; end
        take(k)=take(k)+1;
    end
    while sum(take) > targetN
        eligible = find(take>minPerGroup);
        if isempty(eligible), break; end
        [~,j] = max(take(eligible)-raw(eligible));
        take(eligible(j)) = take(eligible(j))-1;
    end
    idx = zeros(sum(take),1); at = 0;
    for ig = 1:groupCount
        members = find(groupId==ig);
        [~,ord] = sortrows([map.x_v(members), map.t_v(members)],[1 2]);
        members = members(ord);
        q = round(linspace(1,numel(members),take(ig)));
        q = unique(max(1,min(numel(members),q)));
        sel = members(q(:));
        idx(at+(1:numel(sel))) = sel;
        at = at+numel(sel);
    end
    idx = sort(idx(1:at));
    mode = "sensor_turn_stratified";
end
fields = fieldnames(map);
for i = 1:numel(fields)
    value = map.(fields{i});
    if (isnumeric(value) || islogical(value)) && isvector(value) && numel(value) == n
        map.(fields{i}) = value(idx);
    end
end
map.sampling_info = struct('mode', mode, 'requested', targetN, ...
    'used', numel(idx), 'group_count', groupCount);
end

function map = subset_map_uniform(map, targetN)
n = numel(map.t_v);
if n <= targetN
    map.sampling_info = struct('mode', "complete", 'requested', targetN, ...
        'used', n, 'group_count', NaN);
    return;
end
idx = unique(round(linspace(1, n, targetN))).';
names = fieldnames(map);
for k = 1:numel(names)
    value = map.(names{k});
    if (isnumeric(value) || islogical(value) || isstring(value) || iscell(value)) && ...
            numel(value) == n && ~isscalar(value)
        map.(names{k}) = value(idx);
    end
end
map.sampling_info = struct('mode', "uniform", 'requested', targetN, ...
    'used', numel(idx), 'group_count', NaN);
end

function count = count_groups(map)
if isfield(map,'rev_v') && isfield(map,'S_v')
    count = size(unique([map.rev_v(:),map.S_v(:)],'rows'),1);
else
    count = NaN;
end
end

function snrDb = read_snr_db(cfg)
names = {'snrDb','snr_dB','SNR_dB','snr_db'};
snrDb = NaN;
for i = 1:numel(names)
    if isfield(cfg, names{i}) && isscalar(cfg.(names{i}))
        snrDb = cfg.(names{i});
        return;
    end
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
