function trustOptions = make_gap_trust_options(cfg)
%make_gap_trust_options  Collect static-library trust-window settings.

trustOptions = struct();
trustOptions.enable = get_cfg_field(cfg, 'gapTrustEnable', true);
trustOptions.minHalfWidth = get_cfg_field(cfg, 'gapTrustMinHalfWidth', 2.5);
trustOptions.stepHalfWidth = get_cfg_field(cfg, 'gapTrustStepHalfWidth', 0.25);
trustOptions.maxHalfWidth = get_cfg_field(cfg, 'gapTrustMaxHalfWidth', []);
trustOptions.scoreTolerance = get_cfg_field(cfg, 'gapTrustScoreTolerance', 0.08);
trustOptions.derivativeWeight = get_cfg_field(cfg, 'gapTrustDerivativeWeight', 0.20);
trustOptions.sensitivityWeight = get_cfg_field(cfg, 'gapTrustSensitivityWeight', 0.30);
trustOptions.narrowPenaltyWeight = get_cfg_field(cfg, 'gapTrustNarrowPenaltyWeight', 0.03);
trustOptions.minResponseFraction = get_cfg_field(cfg, 'gapTrustMinResponseFraction', 0.12);
trustOptions.minSensitivityFraction = get_cfg_field(cfg, 'gapTrustMinSensitivityFraction', 0.12);
trustOptions.minGridPoints = get_cfg_field(cfg, 'gapTrustMinGridPoints', 40);
trustOptions.selectionMode = "leave_one_gap_static_validation";
end

function val = get_cfg_field(cfg, name, defaultVal)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = cfg.(name);
else
    val = defaultVal;
end
end
