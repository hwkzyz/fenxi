function R = run_inv_log_2_main_method(highMap, templateLib, cfg, staticState)
%run_inv_log_2_main_method  Unified entry for the current default main method.
%
% The project default is Route 30: low-speed gap calibration followed by a
% broad high-speed gap profile, voltage-domain VP candidate generation,
% complete-voltage refinement, and single/dual model-order selection.
% staticState should normally come from the low-speed no-vibration record.
% Historical methods remain available through their original function names.

if nargin < 4
    staticState = [];
end

mainDir = fileparts(mfilename('fullpath'));
localDir = fullfile(mainDir, 'local_func');
if isempty(which('run_low_high_gap_fast_vp_method'))
    addpath(localDir, '-begin');
end

structureMode=lower(string(get_field(cfg,'route30FrequencyStructureMode','general')));
if structureMode=="general"
    R = run_low_high_gap_fast_vp_method(highMap, templateLib, cfg, staticState);
elseif structureMode=="unknown"
    R = run_inv_log_2_unknown_structure(highMap,templateLib,cfg,staticState);
else
    allowed=["single_sync","single_async","dual_sync_sync","dual_sync_async"];
    if ~any(structureMode==allowed)
        error('invlog2:UnknownFrequencyStructureMode',...
            'Unsupported cfg.route30FrequencyStructureMode: %s',structureMode);
    end
    if get_field(cfg,'route30StructuredFastBudgets',true)
        cfg=apply_structured_fast_budgets(cfg,structureMode);
    end
    structuredTic=tic;
    R=run_inv_log_2_structured_main(staticState,highMap,templateLib,cfg,structureMode);
    R.elapsed_s=toc(structuredTic);
end
R.solver_entry = "run_inv_log_2_main_method";
R.solver_version = "route30_main_v2";
if isfield(templateLib, 'gapInterpBasisName')
    R.response_model = string(templateLib.gapInterpBasisName);
elseif isfield(templateLib, 'gapInterpMode')
    R.response_model = string(templateLib.gapInterpMode);
else
    R.response_model = "unknown";
end

function cfg=apply_structured_fast_budgets(cfg,mode)
% Explicit settings supplied by the caller always win.
switch mode
    case "dual_sync_sync"
        cfg.dualSyncGapCount=get_field(cfg,'dualSyncGapCount',15);
        cfg.dualSyncKeepPerGap=get_field(cfg,'dualSyncKeepPerGap',40);
        cfg.dualSyncIteratedReplayCount=get_field(cfg,'dualSyncIteratedReplayCount',400);
        cfg.dualSyncRefineCount=get_field(cfg,'dualSyncRefineCount',12);
    case "dual_sync_async"
        cfg.syncAsyncGapCount=get_field(cfg,'syncAsyncGapCount',11);
        cfg.syncAsyncFrequencyStepHz=get_field(cfg,'syncAsyncFrequencyStepHz',10);
        cfg.syncAsyncKeepPerGap=get_field(cfg,'syncAsyncKeepPerGap',40);
        cfg.syncAsyncIteratedReplayCount=get_field(cfg,'syncAsyncIteratedReplayCount',400);
        cfg.syncAsyncRefineCount=get_field(cfg,'syncAsyncRefineCount',12);
end
end

function v=get_field(s,n,d)
if isstruct(s)&&isfield(s,n)&&~isempty(s.(n)),v=s.(n);else,v=d;end
end
end
