function T = R5_Load_Main07_LowSpeedTemplate(sourceFile)
%R5_LOAD_MAIN07_LOWSPEEDTEMPLATE Load the frozen Main07 low-speed template.

if nargin < 1 || isempty(sourceFile) || ~isfile(sourceFile)
    error('R5:MissingTemplate', 'A valid Main07 low-speed template MAT file is required.');
end
raw = load(sourceFile);
if ~isfield(raw, 'Template')
    error('R5:TemplateSchema', 'Template variable is missing in %s.', sourceFile);
end
T = raw.Template;
T.sourceFile = sourceFile;
if isfield(T, 'Sensor_IDs')
    T.sensor_ids = double(T.Sensor_IDs(:).');
elseif isfield(T, 'SensorBlade')
    T.sensor_ids = unique(arrayfun(@(x) double(x.sensor_id), T.SensorBlade));
else
    T.sensor_ids = [];
end
if isfield(T, 'UnifiedMethodContract')
    T.opr_contract = T.UnifiedMethodContract;
end
end
