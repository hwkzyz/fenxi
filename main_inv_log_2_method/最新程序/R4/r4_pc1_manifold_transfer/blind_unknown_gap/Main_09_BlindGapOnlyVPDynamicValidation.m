function Result=Main_09_BlindGapOnlyVPDynamicValidation(opts)
%MAIN_09_BLINDGAPONLYVPDYNAMICVALIDATION Canonical R4 entry point.
% The historical implementation is archived under main_inv_log_2_method/旧版程序.
if nargin<1, opts=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'local_func'),'-begin');
addpath(fileparts(mfilename('fullpath')),'-begin');
assert_latest_formal_path(root);
Result=Main_10_UnifiedBlindBackendValidation(opts);
end
