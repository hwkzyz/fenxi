function D = R5_Load_Main08_DynamicMap(sourceFile)
%R5_LOAD_MAIN08_DYNAMICMAP Load a Main08 DynamicMap without modifying it.

if nargin < 1 || isempty(sourceFile) || ~isfile(sourceFile)
    error('R5:MissingDynamicMap', 'A valid Main08 DynamicMap MAT file is required.');
end
raw = load(sourceFile);
if ~isfield(raw, 'DynamicMap')
    error('R5:DynamicMapSchema', 'DynamicMap variable is missing in %s.', sourceFile);
end
D = raw.DynamicMap;
D.sourceFile = sourceFile;
if ~isfield(D, 'Window') || isempty(D.Window)
    error('R5:DynamicMapSchema', 'DynamicMap.Window is empty in %s.', sourceFile);
end
end
