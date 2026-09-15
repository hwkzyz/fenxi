function [Template, Diagnostics, outputFile] = Prepare_AdaptiveSG_Template_20251222(forceRebuild)
%PREPARE_ADAPTIVESG_TEMPLATE_20251222 Build the formal adaptive-SG template.
if nargin < 1, forceRebuild = false; end
cfg = Setup_Paths_20251222();
outputFile = cfg.files.lowSpeedTemplate;
if isfile(outputFile) && ~forceRebuild
    S = load(outputFile, 'Template', 'Diagnostics');
    Template = S.Template;
    Diagnostics = S.Diagnostics;
    return;
end
source = load(cfg.files.lowSpeedTemplateSource, 'point_cloud');
base = load(cfg.files.lowSpeedTemplateLegacy, 'Template');
buildCfg = struct('candidateWindowMm', cfg.lowSpeed.candidateWindowMm, ...
    'foldCount', cfg.lowSpeed.groupedFoldCount, ...
    'minBinCount', cfg.lowSpeed.minBinCount);
[Template, Diagnostics] = step07jcore.build_support_aware_adaptive_sg_template( ...
    source.point_cloud, base.Template, buildCfg);
Template.CreatedBy = mfilename;
Template.CreatedOn = datestr(now, 31);
Template.UnifiedMethodContract = cfg.methodFreezeId;
save(outputFile, 'Template', 'Diagnostics', '-v7.3');
writetable(Diagnostics.selection, fullfile(fileparts(outputFile), ...
    'AdaptiveSG_Selection_20251222.csv'));
writetable(Diagnostics.sensor, fullfile(fileparts(outputFile), ...
    'AdaptiveSG_SensorSummary_20251222.csv'));
fprintf('Saved adaptive-SG template: %s\n', outputFile);
end
