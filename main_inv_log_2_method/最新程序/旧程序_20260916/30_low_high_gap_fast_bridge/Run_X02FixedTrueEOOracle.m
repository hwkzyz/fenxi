function [Detail,Summary]=Run_X02FixedTrueEOOracle(snrDb,seeds)
%RUN_X02FIXEDTRUEEOORACLE Conditional continuous-parameter oracle for X02.
% EO=(15,18) is fixed; all calibration, DAQ noise and support rules remain
% identical to the post-freeze support-aware validation.
if nargin<1||isempty(snrDb),snrDb=[25 15 5];end
if nargin<2||isempty(seeds),seeds=1:20;end
Detail=table();Summary=table();
for s=snrDb(:).'
    [D,S]=Run_DualSyncFunnelV1Generalization(s,seeds,'X02_corner_lowamp',true,40,...
        "v2_all_joint1",24,"x02_fixed_true_oracle","support_aware");
    Detail=[Detail;D]; %#ok<AGROW>
    Summary=[Summary;S]; %#ok<AGROW>
end
outDir=fullfile(fileparts(mfilename('fullpath')),'output','x02_fixed_true_oracle');
if ~exist(outDir,'dir'),mkdir(outDir);end
writetable(Detail,fullfile(outDir,'detail.csv'));writetable(Summary,fullfile(outDir,'summary.csv'));
end
