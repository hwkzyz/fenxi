function ctx = bootstrap_inv_log_2_path(callerFile)
%bootstrap_inv_log_2_path  Ensure main_inv_log_2_method is on the MATLAB path.
%
% Subfolder scripts should call this before using any shared helper in the
% main folder. This keeps standalone execution stable even when MATLAB's
% current folder is somewhere else.

if nargin < 1 || isempty(callerFile)
    callerFile = mfilename('fullpath');
end

callerDir = fileparts(callerFile);
mainDir = callerDir;
for ii = 1:8
    if exist(fullfile(mainDir, 'local_func'), 'dir') == 7 && ...
            exist(fullfile(mainDir, 'load_inv_log_2_project_context.m'), 'file') == 2
        break;
    end
    parentDir = fileparts(mainDir);
    if isempty(parentDir) || strcmp(parentDir, mainDir)
        break;
    end
    mainDir = parentDir;
end

addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);
end
