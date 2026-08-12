function cfgOut = define_frequency_range(cfgIn, frequencyRangeHz, frequencyStepHz)
%DEFINE_FREQUENCY_RANGE  Store a physical frequency range and its EO range.
% The conversion is evaluated at the current high-speed rotational rate:
%     EO = f / f_rot,  f_rot = RPM_high / 60.
% This helper only defines the search metadata; it does not silently alter
% the legacy f1Grid/f2Grid used by the production solver.
if nargin < 2 || isempty(frequencyRangeHz), frequencyRangeHz = [100, 1600]; end
if nargin < 3 || isempty(frequencyStepHz), frequencyStepHz = 5; end
frequencyRangeHz = sort(double(frequencyRangeHz(:).'));
if numel(frequencyRangeHz) ~= 2 || any(~isfinite(frequencyRangeHz)) || frequencyRangeHz(1) <= 0
    error('frequencyRangeHz must be a positive two-element range.');
end
if ~isfield(cfgIn, 'RPM_high') || ~isfinite(cfgIn.RPM_high) || cfgIn.RPM_high <= 0
    error('cfg.RPM_high must be a positive scalar.');
end
fRot = cfgIn.RPM_high / 60;
eoRange = frequencyRangeHz / fRot;
eoGrid = (ceil(eoRange(1)):floor(eoRange(2))).';
fGrid = (frequencyRangeHz(1):frequencyStepHz:frequencyRangeHz(2)).';
if fGrid(end) < frequencyRangeHz(2), fGrid(end+1) = frequencyRangeHz(2); end
cfgOut = cfgIn;
cfgOut.frequencyRangeHz = frequencyRangeHz;
cfgOut.frequencyStepHz = frequencyStepHz;
cfgOut.rotationalFrequencyHz = fRot;
cfgOut.frequencyRangeEO = eoRange;
cfgOut.synchronousEOGrid = eoGrid;
cfgOut.frequencyGridFromRangeHz = fGrid;
end
