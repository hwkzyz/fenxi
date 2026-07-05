function cfg = BTTProjectConfig_20250527()
%BTTPROJECTCONFIG_20250527 Project-level configuration for the 20250527 BTT route.
%
% This wrapper keeps the Step scripts semantically clear: cfg is for stable
% project/data facts, while each Step keeps its own local settings at the top
% of the script. The underlying data facts remain in BTTDataConfig_20250527
% for compatibility with existing scripts.

cfg = BTTDataConfig_20250527();
end
