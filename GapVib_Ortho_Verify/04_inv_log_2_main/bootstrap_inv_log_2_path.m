function ctx = bootstrap_inv_log_2_path(callerFile)
%bootstrap_inv_log_2_path  Ensure 04_inv_log_2_main is on the MATLAB path.
%
% Subfolder scripts should call this before using any shared helper in the
% main folder. This keeps standalone execution stable even when MATLAB's
% current folder is somewhere else.

if nargin < 1 || isempty(callerFile)
    callerFile = mfilename('fullpath');
end

callerDir = fileparts(callerFile);
mainDir = fileparts(callerDir);
addpath(mainDir, '-begin');
ctx = load_inv_log_2_project_context(mainDir);
end
