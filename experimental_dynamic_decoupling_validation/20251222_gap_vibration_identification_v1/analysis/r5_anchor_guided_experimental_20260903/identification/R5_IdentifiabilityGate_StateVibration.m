function gate = R5_IdentifiabilityGate_StateVibration(profile, candidates, cfg)
%R5_IDENTIFIABILITYGATE_STATEVIBRATION Reject state/vibration confounding.
P = profile.stateProfiles;
valid = strcmp({P.status},'profile_fitted');
Pv = P(valid);
gate = struct('status','rejected','profileMarginMv2',NaN,'rho_sx',NaN, ...
    'jacobianCond',Inf,'amplitudeAtBound',false,'dxAtBound',false,'reasons',{{}});
if isempty(Pv), gate.reasons={'no_profile_fit'}; return; end
J = sort([Pv.J]);
if numel(J)>1, gate.profileMarginMv2=J(2)-J(1); else, gate.profileMarginMv2=Inf; end
s = profile.selected.stateValue;
states = [candidates.stateValue]; [~,k] = min(abs(states-s));
k0=max(1,k-1); k1=min(numel(candidates),k+1);
x = candidates(k).x_mm(:); dFdx=gradient(candidates(k).waveform_mv(:),x);
dFds=(candidates(k1).waveform_mv(:)-candidates(k0).waveform_mv(:)) / max(states(k1)-states(k0),eps);
mask=isfinite(dFdx)&isfinite(dFds); R=corrcoef(dFdx(mask),dFds(mask));
if numel(R)==4, gate.rho_sx=R(1,2); end
G=[dFdx(mask),dFds(mask)]; gate.jacobianCond=cond(G.'*G);
if gate.profileMarginMv2 < cfg.gates.minProfileMarginMv2, gate.reasons{end+1}='profile_margin'; end
if abs(gate.rho_sx) > cfg.gates.maxAbsRhoSx, gate.reasons{end+1}='state_vibration_correlation'; end
if gate.jacobianCond > cfg.gates.maxJacobianCondition, gate.reasons{end+1}='jacobian_condition'; end
tolA = max(1e-5, 0.01 * diff(cfg.profile.amplitudeBoundsMm));
tolX = max(1e-5, 0.01 * diff(cfg.profile.dxBoundsMm));
gate.amplitudeAtBound = abs(profile.selected.amplitudeMm-cfg.profile.amplitudeBoundsMm(1)) <= tolA || ...
    abs(profile.selected.amplitudeMm-cfg.profile.amplitudeBoundsMm(2)) <= tolA;
gate.dxAtBound = abs(profile.selected.dxMm-cfg.profile.dxBoundsMm(1)) <= tolX || ...
    abs(profile.selected.dxMm-cfg.profile.dxBoundsMm(2)) <= tolX;
if gate.amplitudeAtBound, gate.reasons{end+1}='amplitude_bound'; end
if gate.dxAtBound, gate.reasons{end+1}='dx_bound'; end
% This is a per-window mathematical gate only.  Deployment additionally
% requires the external V1/V2/V5 acceptance criteria and sensor agreement.
if isempty(gate.reasons), gate.status='identifiable_window_only'; else, gate.status='diagnostic_only'; end
end
