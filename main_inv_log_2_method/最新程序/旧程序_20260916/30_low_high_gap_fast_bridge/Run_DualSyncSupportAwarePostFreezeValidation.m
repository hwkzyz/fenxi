function [Detail,Summary]=Run_DualSyncSupportAwarePostFreezeValidation(snrDb,seeds,conditionIds)
%RUN_DUALSYNCSUPPORTAWAREPOSTFREEZEVALIDATION Frozen V2 support validation.
% The physical envelope is declared in Run_DualSyncFunnelV1Generalization
% from paper-wide simulation ranges, independently of X01 outcomes.
if nargin<1||isempty(snrDb),snrDb=[25 15 5];end
if nargin<2||isempty(seeds),seeds=1:20;end
if nargin<3,conditionIds=[];end
Detail=table();Summary=table();
for level=snrDb(:).'
    [D,S]=Run_DualSyncFunnelV1Generalization(level,seeds,conditionIds,false,40,...
        "v2_all_joint1",24,"postfreeze_unconstrained","support_aware");
    Detail=[Detail;D]; %#ok<AGROW>
    Summary=[Summary;S]; %#ok<AGROW>
end
outDir=fullfile(fileparts(mfilename('fullpath')),'output','dual_sync_supportaware_postfreeze');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Detail,fullfile(outDir,'detail.csv'));writetable(Summary,fullfile(outDir,'summary.csv'));
disp(Summary(:,intersect({'condition_id','nominal_snr_db','n','eo_correct_rate',...
    'success_rate','support_runtime_valid_rate'},Summary.Properties.VariableNames,'stable')));
end
