function cfg = make_identifiability_analysis_config()
%MAKE_IDENTIFIABILITY_ANALYSIS_CONFIG Fixed settings for mechanism analysis.
cfg = struct();
cfg.version = "identifiability_mechanism_analysis_v1";
cfg.identity_tolerance = 1e-10;
cfg.analysis_scope = "P0-P3 minimum falsifiable loop";
cfg.modify_frozen_solver = false;
cfg.tune_funnel_or_template = false;
end

