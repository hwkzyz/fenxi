function R = run_two_vp_gap_vib_joint(highMap, templateLib, cfg, staticState)
%run_two_vp_gap_vib_joint  Fast gap init -> VP-aided gap correction -> VP vibration fit -> joint refine.
%
% The first VP pass estimates the vibration-shaped nuisance subspace used
% to correct the gap. The second VP pass identifies vibration parameters
% under the corrected gap before final full-model nonlinear refinement.

t = highMap.t_v(:);
x = highMap.x_v(:);
V = highMap.V_a(:);
[~, ~, gidx] = unique([highMap.rev_v(:), highMap.S_v(:)], 'rows');
ng = max(gidx);

if nargin < 4 || isempty(staticState) || ~isfield(staticState, 'gHat')
    staticState = estimate_highspeed_static_gap_raw(highMap, templateLib, cfg);
end

timerGapCorr = tic;
cfgGap = cfg;
cfgGap.gapInitFreqMode = "vp_fixed_gap";
gapStateProjected = estimate_gap_init_vib_basis_projected(highMap, templateLib, cfgGap, staticState);
timeGapCorr = toc(timerGapCorr);

muList = get_cfg_field(cfg, 'twoVpJointMuList', [0, 0.15, 0.75]);
timerBranch = tic;
branchProjected = evaluate_gap_branch("vib_projected_gap", gapStateProjected, t, V, x, templateLib, cfg, gidx, ng, muList);
timeBranch = toc(timerBranch);

bestBranch = branchProjected;
bestRmse = bestBranch.rmse;
bestJoint = bestBranch.joint;
fit1 = bestBranch.fit;
gapState = bestBranch.gapState;

VFit = eval_gap_template(templateLib, bestJoint.g, x - bestJoint.dx - fit_u(bestJoint.p, t));

R.method = "two_vp_gap_vib_joint";
R.g_used = bestJoint.g;
R.dx_used = bestJoint.dx;
R.fit = bestJoint;
R.f_id = bestJoint.f;
R.A_id = bestJoint.A;
R.phi_id = bestJoint.phi;
R.VFit = VFit;
R.rmse = bestRmse;
R.eta_g = eta_gap(templateLib, bestJoint.g);
R.deltaJ = fit1.deltaJ;
R.rhoJ = fit1.rhoJ;
R.staticState = staticState;
R.gapState = gapState;
R.used_gap_correction = bestBranch.name == "vib_projected_gap";
R.selected_branch = bestBranch.name;
if isfield(gapState, 'freq_meta') && isfield(gapState.freq_meta, 'rhoJ')
    R.first_vp_rhoJ = gapState.freq_meta.rhoJ;
else
    R.first_vp_rhoJ = NaN;
end
R.vpFit = fit1;
R.time_gap_correction_s = timeGapCorr;
R.time_vp_fit_s = bestBranch.time_vp_s;
R.time_joint_s = bestBranch.time_joint_s;
R.time_branch_competition_s = timeBranch;
R.branch_table = struct2table(rmfield(branchProjected, {'joint','fit','gapState'}));
end

function branch = evaluate_gap_branch(name, gapState, t, V, x, templateLib, cfg, gidx, ng, muList)
timerVp = tic;
fit0 = fit_vp_main(t, V, x, templateLib, gapState.gHat, cfg);
topK = get_cfg_field(cfg, 'twoVpTopKFreqCandidates', 5);
candidateList = fit0.candidates(1:min(topK, numel(fit0.candidates)));
refinedList = repmat(struct('fit', [], 'sse', inf, 'candidate_index', 0), numel(candidateList), 1);
for ic = 1:numel(candidateList)
    fitCand = candidateList(ic);
    fitCand.J1 = fit0.J1;
    fitCand.J2 = fit0.J2;
    fitCand.deltaJ = fit0.deltaJ;
    fitCand.rhoJ = fit0.rhoJ;
    fitTry = refine_vp_main_fit(fitCand, t, V, x, templateLib, gapState.gHat, cfg);
    resTry = V - eval_gap_template(templateLib, gapState.gHat, x - fitTry.delta_sample);
    refinedList(ic).fit = fitTry;
    refinedList(ic).sse = dot(resTry, resTry);
    refinedList(ic).candidate_index = ic;
end
timeVp = toc(timerVp);

[~, refineOrder] = sort([refinedList.sse], 'ascend');
jointKeep = min(get_cfg_field(cfg, 'twoVpJointCandidateKeep', 2), numel(refineOrder));

bestJoint = [];
bestRmse = inf;
bestFit = [];
timerJoint = tic;
for ik = 1:jointKeep
    refIdx = refineOrder(ik);
    fitTry = refinedList(refIdx).fit;
    for imu = 1:numel(muList)
        jointTry = refine_projected_joint(t, V, x, templateLib, cfg, gidx, ng, ...
            gapState.gHat, gapState.dx0, fitTry.p, muList(imu));
        VTry = eval_gap_template(templateLib, jointTry.g, ...
            x - jointTry.dx - fit_u(jointTry.p, t));
        rmseTry = sqrt(mean((V - VTry).^2));
        jointTry.mu = muList(imu);
        jointTry.wave_rmse = rmseTry;
        jointTry.candidate_index = refinedList(refIdx).candidate_index;
        if rmseTry < bestRmse
            bestRmse = rmseTry;
            bestJoint = jointTry;
            bestFit = fitTry;
        end
    end
end
timeJoint = toc(timerJoint);

branch = struct();
branch.name = string(name);
branch.rmse = bestRmse;
branch.gap_start = gapState.gHat;
branch.gap_final = bestJoint.g;
branch.f1 = bestJoint.f(1);
branch.f2 = bestJoint.f(2);
branch.deltaJ = bestFit.deltaJ;
branch.rhoJ = bestFit.rhoJ;
branch.candidate_index = bestJoint.candidate_index;
branch.candidate_count = numel(candidateList);
branch.joint_candidate_keep = jointKeep;
branch.time_vp_s = timeVp;
branch.time_joint_s = timeJoint;
branch.joint = bestJoint;
branch.fit = bestFit;
branch.gapState = gapState;
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
