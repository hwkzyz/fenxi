function filePath = resolve_project_data_file(rootDir, savedPath)
%resolve_project_data_file  Resolve saved data paths after workspace cleanup.

if exist(savedPath, 'file')
    filePath = savedPath;
    return;
end

[~, fileName, ext] = fileparts(savedPath);
targetName = [fileName, ext];
candidates = {
    fullfile(rootDir, 'data', '间隙的影响', targetName)
    fullfile(rootDir, '间隙的影响', targetName)
    fullfile(rootDir, 'data', targetName)
    };

for ii = 1:numel(candidates)
    if exist(candidates{ii}, 'file')
        filePath = candidates{ii};
        return;
    end
end

error('Could not locate data file: %s', savedPath);
end
