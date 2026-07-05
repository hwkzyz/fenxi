%% Step07_AuditNumberingAndCorrespondence_20251222
% Lightweight audit of NewFlow region, numbering, template, and ID outputs.

clear; clc;
P = NewFlow_Config_20251222();
ensure_parent_dir_local(P.files.auditSummary);

required = {
    P.files.regionSelection, 'RegionSelection'
    P.files.highSpeedNumbering, 'HighSpeedNumbering'
    P.files.lowSpeedTemplateLibrary, 'LowSpeedTemplateLibrary'
    P.files.identificationResult, 'IdentificationResult'
    P.files.identificationSummary, 'IdentificationSummary'
    };

rows = repmat(struct('Artifact', '', 'Path', '', 'Exists', false, 'Bytes', NaN), ...
    size(required, 1), 1);
for i = 1:size(required, 1)
    pathText = required{i, 1};
    info = dir(pathText);
    rows(i).Artifact = required{i, 2};
    rows(i).Path = pathText;
    rows(i).Exists = ~isempty(info);
    if ~isempty(info)
        rows(i).Bytes = info(1).bytes;
    end
end

AuditSummary = struct2table(rows);
writetable(AuditSummary, P.files.auditSummary);

fprintf('\n=== Step07: NewFlow audit ===\n');
disp(AuditSummary(:, {'Artifact','Exists','Bytes'}));
fprintf('Saved: %s\n', P.files.auditSummary);

function ensure_parent_dir_local(pathText)
folder = fileparts(pathText);
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
end
