function [Detail,Summary]=Run_DualSyncFunnelV2Generalization(snrDb,seeds,conditionIds,secondJointK,supportMode)
%RUN_DUALSYNCFUNNELV2GENERALIZATION Formal dual-synchronous Funnel V2 entry.
if nargin<1||isempty(snrDb),snrDb=15;end
if nargin<2||isempty(seeds),seeds=1:10;end
if nargin<3,conditionIds=[];end
if nargin<4||isempty(secondJointK),secondJointK=24;end
if nargin<5||isempty(supportMode),supportMode="legacy";end
[Detail,Summary]=Run_DualSyncFunnelV1Generalization(...
    snrDb,seeds,conditionIds,false,40,"v2_all_joint1",secondJointK,"",supportMode);
end
