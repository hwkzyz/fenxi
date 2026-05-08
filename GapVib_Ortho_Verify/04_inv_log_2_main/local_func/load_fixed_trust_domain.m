function trustInfo = load_fixed_trust_domain(thisDir)
%load_fixed_trust_domain  Load a previously calibrated trust domain.

trustPath = get_trust_result_path(thisDir);
if ~exist(trustPath, 'file')
    error(['Fixed trust-domain file is missing:\n%s\n', ...
        'Run 01_trust_domain_calibration/Step00_Calibrate_TrustDomain.m first.'], trustPath);
end
S = load(trustPath, 'trustInfo');
trustInfo = S.trustInfo;
end
