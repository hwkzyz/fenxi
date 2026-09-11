function [file, candidates] = R5_Find_Main08_DynamicMap(cfg, explicitFile)
%R5_FIND_MAIN08_DYNAMICMAP Find an existing Main08 output without rebuilding it.
if nargin < 2, explicitFile = ''; end
if ~isempty(explicitFile)
    if ~isfile(explicitFile), error('R5:DynamicMapMissing', 'Missing DynamicMap: %s', explicitFile); end
    file = explicitFile; candidates = {explicitFile}; return
end
root = cfg.productionRoot;
d = dir(fullfile(root, '**', 'Step06_DynamicMap_*.mat'));
candidates = arrayfun(@(q) fullfile(q.folder, q.name), d, 'UniformOutput', false);
if isempty(candidates)
    file = '';
    return
end
if numel(candidates) > 1
    error('R5:DynamicMapAmbiguous', 'Found %d DynamicMap files; pass an explicit file to R5_Load_Main08_DynamicMap.', numel(candidates));
end
file = candidates{1};
end
