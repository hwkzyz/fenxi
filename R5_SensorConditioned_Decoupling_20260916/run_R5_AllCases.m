function out = run_R5_AllCases()
%RUN_R5_ALLCASES Run both packaged formal R5 cases.
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir,'-begin');
results = struct();
results.case20250527 = run_R5_20250527();
results.case20251222 = run_R5_20251222();
out = results;
end
