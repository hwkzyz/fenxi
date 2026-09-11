function cfg = BTTProjectConfig_20241106()
%BTTPROJECTCONFIG_20241106 Project-level wrapper for the 20241106 route.
%
% Keep this wrapper so later method scripts can ask for project config while
% the stable data facts remain centralized in BTTDataConfig_20241106.

cfg = BTTDataConfig_20241106();
end
