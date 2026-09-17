function Result=Run_StructuredFourModeRepeatRegression(seedOffsets)
%RUN_STRUCTUREDFOURMODEREPEATREGRESSION Repeat the 5 dB acceptance across noise seeds.
if nargin<1,seedOffsets=[0 1000 2000];end
allAudit=table();
for i=1:numel(seedOffsets)
    q=Run_StructuredFourModeRegression(seedOffsets(i));
    A=q.Audit;A.seed_offset(:)=seedOffsets(i);A.minimum_revolution_gate(:)=q.minimumRevolutionGatePassed;
    allAudit=[allAudit;A]; %#ok<AGROW>
end
root=fileparts(mfilename('fullpath'));out=fullfile(root,'output','structured_four_mode_repeat');
if ~exist(out,'dir'),mkdir(out);end
writetable(allAudit,fullfile(out,'structured_four_mode_repeat.csv'));
save(fullfile(out,'structured_four_mode_repeat.mat'),'allAudit','seedOffsets');
summary=groupsummary(allAudit,'mode',{'sum','mean','max'},...
    {'success','identification_confident','max_frequency_error','elapsed_s'});
writetable(summary,fullfile(out,'structured_four_mode_repeat_summary.csv'));
disp(summary);
if ~all(allAudit.success),error('Repeated four-mode structured regression failed.');end
Result=struct('Audit',allAudit,'Summary',summary,'outputDir',out);
end
