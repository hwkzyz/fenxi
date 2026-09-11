function [eta, info] = load_joint_static_eta_preview(foundationDir, bladeId, sensorIds, varargin)
%LOAD_JOINT_STATIC_ETA_PREVIEW Load the fixed joint static eta used by Step05.

opts = struct('StaticEtaFile', '');
if mod(numel(varargin), 2) ~= 0
    error('Optional arguments must be name/value pairs.');
end
for i = 1:2:numel(varargin)
    name = char(varargin{i});
    value = varargin{i + 1};
    if strcmpi(name, 'StaticEtaFile') || strcmpi(name, 'EtaFile')
        opts.StaticEtaFile = char(value);
    else
        error('Unsupported option: %s', name);
    end
end

sensorIds = sensorIds(:).';
bladeId = double(bladeId);
sensorTag = ['S', sprintf('%d', sensorIds)];

etaFile = strtrim(opts.StaticEtaFile);
if isempty(etaFile)
    etaFile = fullfile(foundationDir, 'output', ...
        'step05_joint_static_eta_preview', ...
        sprintf('JointStaticEtaPreview_B%d_%s.mat', bladeId, sensorTag));
end

if ~isfile(etaFile)
    error('Missing joint static eta preview: %s', etaFile);
end

loaded = load(etaFile, 'Summary');
if ~isfield(loaded, 'Summary')
    error('Joint static eta file must contain Summary: %s', etaFile);
end
Summary = loaded.Summary;

if isfield(Summary, 'TargetBlade') && double(Summary.TargetBlade) ~= bladeId
    error('Joint static eta target blade mismatch: file B%d, requested B%d.', ...
        double(Summary.TargetBlade), bladeId);
end
if isfield(Summary, 'SensorIDs') && ...
        ~isequal(double(Summary.SensorIDs(:).'), sensorIds)
    error('Joint static eta sensor mismatch: file %s, requested %s.', ...
        mat2str(double(Summary.SensorIDs(:).')), mat2str(sensorIds));
end
if ~isfield(Summary, 'Consensus')
    error('Joint static eta file missing Summary.Consensus: %s', etaFile);
end

consensus = Summary.Consensus;
if isfield(consensus, 'eta_median_mm')
    eta = double(consensus.eta_median_mm(:).');
    etaField = 'eta_median_mm';
elseif isfield(consensus, 'eta_plan_mm')
    eta = double(consensus.eta_plan_mm(:).');
    etaField = 'eta_plan_mm';
else
    error(['Joint static eta file must contain Summary.Consensus.' ...
        'eta_median_mm or eta_plan_mm: %s'], etaFile);
end

if numel(eta) ~= numel(sensorIds) || any(~isfinite(eta))
    error('Invalid joint static eta vector in %s.', etaFile);
end
eta(1) = 0;

planIqr = nan(size(eta));
if isfield(consensus, 'eta_plan_iqr_mm')
    planIqr = double(consensus.eta_plan_iqr_mm(:).');
end
if numel(planIqr) ~= numel(eta)
    planIqr = nan(size(eta));
end

info = struct();
info.file = etaFile;
info.sourceField = etaField;
info.targetBlade = bladeId;
info.sensorIds = sensorIds;
info.sensorTag = sensorTag;
info.etaMm = eta;
info.etaPlanIqrMm = planIqr;
info.table = table(sensorIds(:), eta(:), planIqr(:), ...
    repmat(string(etaFile), numel(sensorIds), 1), ...
    'VariableNames', {'sensor_id', 'eta_prior_mm', ...
    'eta_plan_iqr_mm', 'joint_static_eta_file'});
end
