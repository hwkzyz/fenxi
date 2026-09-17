function [pairs, info] = generate_fast_dual_frequency_pairs(map, lib, cfg, lowState, externalAnchors)
%GENERATE_FAST_DUAL_FREQUENCY_PAIRS Conditional batch dual-support screen.
%
% The screen uses only measured time, sensor position, voltage, the forward
% library, and the low-speed gap estimate.  Simulation truth is deliberately
% absent.  Data-driven single-frequency supports anchor batched scans for the
% missing component over a few gap anchors.  This is O(Nanchor*Ngap*Nf), not
% O(Ngap*Nf^2). Low-confidence screens request the full-grid fallback.

t0 = tic;
if nargin<5,externalAnchors=[];end
t = map.t_v(:); x = map.x_v(:); V = map.V_a(:);
fGrid = unique(sort(cfg.f1Grid(:))).';
minSep = get_cfg(cfg, 'route30FastCandidateMinSeparationHz', 20);
nGrid = numel(fGrid);
fullPairCount = nnz(triu(abs(fGrid(:)-fGrid(:).') >= minSep, 1));

gHalf = get_cfg(cfg, 'mainVpGapHalfWidth', .70);
gLo = max(.05, lowState.gHat-gHalf);
gHi = lowState.gHat+gHalf;
if isfield(lib,'gapTrain') && ~isempty(lib.gapTrain)
    margin = get_cfg(cfg,'rawGapSearchMargin',0);
    gLo = max(gLo,min(lib.gapTrain)-margin);
    gHi = min(gHi,max(lib.gapTrain)+margin);
end
nAnchor = max(1,round(get_cfg(cfg,'route30FastCandidateGapAnchors',5)));
gapAnchors = unique(linspace(gLo,gHi,nAnchor));
keepPerGap = max(10,round(get_cfg(cfg,'route30FastCandidatePairsPerGap',16)));

anchorCount=max(1,round(get_cfg(cfg,'route30FastCandidateSingleAnchorCount',8)));
anchorDiversity=get_cfg(cfg,'route30FastCandidateSingleAnchorDiversityHz',10);
externalAnchors=externalAnchors(isfinite(externalAnchors));
singleAnchors=[];
for f=externalAnchors(:).'
    if isempty(singleAnchors)||all(abs(singleAnchors-f)>=anchorDiversity)
        singleAnchors(end+1)=f; %#ok<AGROW>
        if numel(singleAnchors)>=anchorCount,break;end
    end
end

probe = cfg;
probe.f1Grid = fGrid;
probe.vpUniqueUnorderedPairs = false;
probe.vpUseStagedGrid = false;
probe.numVarproCandidates = keepPerGap;
probe.vpCandidateDiversityHz = get_cfg(cfg,'route30FastCandidateDiversityHz',5);
probe.vpCandidateFrequencyPairsHz = [];

pairsRaw = zeros(0,2); scoreRaw=zeros(0,1);
anchorBest = nan(numel(gapAnchors),2);
explained = nan(numel(gapAnchors),1); weakRatio = nan(numel(gapAnchors),1);
anchorTime = zeros(numel(gapAnchors),1);
for ig = 1:numel(gapAnchors)
    ta = tic;
    F0 = eval_gap_template(lib,gapAnchors(ig),x);
    y = V-F0;
    y2=max(dot(y,y),eps);bestScore=Inf;bestFit=[];
    for ia=1:numel(singleAnchors)
        probe.f2Grid=singleAnchors(ia);
        fit = fit_vp_main(t,V,x,lib,gapAnchors(ig),probe);
        q = vertcat(fit.candidates.f);
        s = [fit.candidates.linear_sse].'/y2;
        pairsRaw=[pairsRaw;sort(q,2)];scoreRaw=[scoreRaw;s]; %#ok<AGROW>
        if fit.J1<bestScore,bestScore=fit.J1;bestFit=fit;end
    end
    anchorTime(ig) = toc(ta);
    if ~isempty(bestFit)
        anchorBest(ig,:) = sort(bestFit.candidates(1).f);
        explained(ig) = max(0,1-bestFit.J1/y2);
        a = bestFit.candidates(1).A;
        weakRatio(ig) = min(a)/max(max(a),eps);
    end
end
[~,ord]=sort(scoreRaw,'ascend');pairsRaw=pairsRaw(ord,:);
pairs = unique(round(sort(pairsRaw,2)*1e9)/1e9,'rows','stable');
maxPairs = round(get_cfg(cfg,'route30FastCandidateMaxPairs',640));
if size(pairs,1)>maxPairs, pairs=pairs(1:maxPairs,:); end

minExplained = get_cfg(cfg,'route30FastCandidateMinExplainedFraction',.015);
minWeakRatio = get_cfg(cfg,'route30FastCandidateMinWeakRatio',.04);
minPairCount=round(get_cfg(cfg,'route30FastCandidateMinPairCount',80));
confident = ~isempty(singleAnchors) && size(pairs,1)>=minPairCount && any(explained>=minExplained) && ...
    any(weakRatio>=minWeakRatio) && all(isfinite(pairs),'all');
reason = "accepted";
if isempty(singleAnchors),reason="no_single_frequency_support";
elseif size(pairs,1)<minPairCount, reason="too_few_pairs";
elseif ~any(explained>=minExplained), reason="low_explained_fraction";
elseif ~any(weakRatio>=minWeakRatio), reason="weak_second_support";
elseif ~all(isfinite(pairs),'all'), reason="nonfinite_candidate";
end

info = struct('method',"derivative_weighted_conditional_batch_topk",...
    'confident',confident,'fallback_to_full',~confident,'reason',reason,...
    'frequency_grid_count',nGrid,'full_unordered_pair_count',fullPairCount,...
    'gap_anchor_count',numel(gapAnchors),'gap_anchors_mm',gapAnchors,...
    'single_anchor_frequencies_hz',singleAnchors,...
    'pairs_per_anchor',keepPerGap,'candidate_pair_count',size(pairs,1),...
    'candidate_pairs_hz',pairs,'anchor_best_pairs_hz',anchorBest,...
    'anchor_explained_fraction',explained,'anchor_weak_ratio',weakRatio,...
    'anchor_time_s',anchorTime,'elapsed_s',toc(t0),...
    'uses_truth',false,'library_boundary_touched',...
    any(abs(gapAnchors-min(lib.gapTrain))<1e-12 | ...
        abs(gapAnchors-max(lib.gapTrain))<1e-12));
end

function value=get_cfg(s,name,defaultValue)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),value=s.(name);
else,value=defaultValue;end
end
