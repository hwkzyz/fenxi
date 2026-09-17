function cfg=configure_structured_identification(cfg,profile)
%CONFIGURE_STRUCTURED_IDENTIFICATION Apply verified structured-route settings.
if nargin<2||isempty(profile),profile="experiment_current";end
profile=lower(string(profile));
cfg.structuredMinRevolutions=8;
cfg.structuredRequireMinRevolutions=true;
cfg.structuredConfidenceRevolutions=16;
cfg.structuredMaximumRevolutions=32;
cfg.structuredComponentAmplitudeFloorMm=.075;
cfg.dualSyncKeepPerGap=max(100,get_field(cfg,'dualSyncKeepPerGap',0));
cfg.dualSyncIteratedReplayCount=max(2100,get_field(cfg,'dualSyncIteratedReplayCount',0));
cfg.dualSyncRefineCount=max(50,get_field(cfg,'dualSyncRefineCount',0));
cfg.dualSyncSearchStrategy=get_field(cfg,'dualSyncSearchStrategy',"funnel_v1");
cfg.funnelAllPairsFirstJointStep=get_field(cfg,'funnelAllPairsFirstJointStep',true);
cfg.funnelSecondJointPairCount=get_field(cfg,'funnelSecondJointPairCount',24);
cfg.funnelFullRefineCount=get_field(cfg,'funnelFullRefineCount',3);
cfg.syncAsyncKeepPerGap=max(100,get_field(cfg,'syncAsyncKeepPerGap',0));
cfg.syncAsyncIteratedReplayCount=max(2100,get_field(cfg,'syncAsyncIteratedReplayCount',0));
cfg.syncAsyncRefineCount=max(50,get_field(cfg,'syncAsyncRefineCount',0));
switch profile
    case "experiment_current"
        cfg.structuredLayoutProfile="measured_hardware";
    case "simulation_robust5"
        cfg.alpha_k=[0 .973466 2.42106 2.90560 4.06032];
        cfg.structuredLayoutProfile="eo_subspace_optimized_5_probe";
    otherwise
        error('invlog2:UnknownStructuredProfile','Unknown structured profile: %s',profile);
end
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
