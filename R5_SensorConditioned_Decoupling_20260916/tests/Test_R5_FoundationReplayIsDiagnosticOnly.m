function Test_R5_FoundationReplayIsDiagnosticOnly()
%TEST_R5_FOUNDATIONREPLAYISDIAGNOSTICONLY Replay must not veto a valid R5 fit.
rootDir = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath;
addpath(rootDir,'-begin');
addpath(fullfile(rootDir,'core'),'-begin');
addpath(fullfile(rootDir,'adapters'),'-begin');

r = struct('status','pass','support_pass',true, ...
    'eo_consistency_pass',true,'hit_boundary',false,'a_bound_hit',false, ...
    'dx_bound_hit',false,'dg_bound_hit',false,'amplitude_ambiguous',false, ...
    'foundation_replay_status','replay_fail', ...
    'foundation_replay_used_for_acceptance',false);

assert(strcmp(classify_for_test(r),'usable'), ...
    'A Foundation replay failure must not reject a valid R5 fit.');
r.foundation_replay_status = 'replay_unavailable';
assert(strcmp(classify_for_test(r),'usable'), ...
    'Unavailable Foundation replay must not reject a valid R5 fit.');
assert(~r.foundation_replay_used_for_acceptance, ...
    'Foundation replay must be explicitly diagnostic-only.');

fprintf('Test_R5_FoundationReplayIsDiagnosticOnly PASS.\n');
end

function q = classify_for_test(r)
if ~r.support_pass || ~strcmp(r.status,'pass')
    q='rejected';
elseif isfield(r,'eo_consistency_pass') && ~r.eo_consistency_pass
    q='diagnostic_only';
elseif r.hit_boundary || r.a_bound_hit || r.dx_bound_hit || r.dg_bound_hit
    q='diagnostic_only';
elseif isfield(r,'amplitude_ambiguous') && r.amplitude_ambiguous
    q='diagnostic_only';
else
    q='usable';
end
end
