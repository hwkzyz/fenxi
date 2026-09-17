function M = compute_phase_diversity_metrics(phase,weight)
%COMPUTE_PHASE_DIVERSITY_METRICS Harmonic phase-diversity reference metrics.
phase = phase(:); weight = weight(:);
if numel(phase) ~= numel(weight)
    error('identifiability:SizeMismatch','phase and weight must have equal length.');
end
valid = isfinite(phase) & isfinite(weight) & weight >= 0;
phase = phase(valid); weight = weight(valid); W0 = sum(weight);
if W0 <= 0
    error('identifiability:InvalidWeight','At least one positive weight is required.');
end
r2 = abs(sum(weight.*exp(1i*2*phase)))/W0;
r2 = min(max(real(r2),0),1);
lambdaMin = 0.5*W0*(1-r2); lambdaMax = 0.5*W0*(1+r2);
M = struct('sample_count',numel(phase),'weight_sum',W0,'r2',r2, ...
    'lambda_min_gram',lambdaMin,'lambda_max_gram',lambdaMax, ...
    'determinant_gram',0.25*W0^2*(1-r2^2), ...
    'condition_gram',lambdaMax/max(lambdaMin,eps), ...
    'condition_design',sqrt(lambdaMax/max(lambdaMin,eps)));
end

