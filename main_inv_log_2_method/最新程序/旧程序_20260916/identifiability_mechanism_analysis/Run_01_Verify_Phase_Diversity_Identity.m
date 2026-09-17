function Result = Run_01_Verify_Phase_Diversity_Identity()
%RUN_01_VERIFY_PHASE_DIVERSITY_IDENTITY Validate the analytic r2 identity.
analysisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(analysisDir,'config'),'-begin');
addpath(fullfile(analysisDir,'local_func'),'-begin');
cfg = make_identifiability_analysis_config();
phase = linspace(-0.9,1.1,401).';
weight = 0.2 + exp(-0.5*(phase/0.42).^2);
M = compute_phase_diversity_metrics(phase,weight);
X = sqrt(weight).*[sin(phase),cos(phase)];
lambdaNumeric = sort(eig(X.'*X),'ascend').';
lambdaAnalytic = [M.lambda_min_gram,M.lambda_max_gram];
absoluteError = max(abs(lambdaNumeric-lambdaAnalytic));
if absoluteError > cfg.identity_tolerance
    error('identifiability:PhaseIdentityFailed', ...
        'Identity error %.3g exceeds tolerance %.3g.',absoluteError,cfg.identity_tolerance);
end
Result = struct('metrics',M,'lambda_numeric',lambdaNumeric, ...
    'lambda_analytic',lambdaAnalytic,'max_absolute_error',absoluteError, ...
    'identity_tolerance',cfg.identity_tolerance);
outDir = fullfile(analysisDir,'output','01_phase_diversity_identity');
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir,'phase_diversity_identity.mat'),'Result');
writetable(struct2table(Result.metrics),fullfile(outDir,'phase_diversity_metrics.csv'));
fprintf('Phase-diversity identity passed; max eigenvalue error = %.3g.\n',absoluteError);
end

