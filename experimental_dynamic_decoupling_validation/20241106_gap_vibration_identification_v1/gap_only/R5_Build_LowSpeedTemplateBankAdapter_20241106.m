function outputFile = R5_Build_LowSpeedTemplateBankAdapter_20241106(sourceFile, outputFile)
%R5_BUILD_LOWSPEEDTEMPLATEBANKADAPTER_20241106 Adapt the legacy bank to R5.
% The bank stores all blade/sensor templates in mV.  R5's template contract
% stores baseline-removed values in V, so conversion is explicit here.
cfg = Config_20241106();
if nargin < 1 || isempty(sourceFile)
    sourceFile = fullfile(cfg.paths.calibrationInputs,'gap', ...
        'LowSpeedTemplateBank_20241106_B1toB6_S2357.mat');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(cfg.paths.results,'r5_low_speed_template_bank_adapter.mat');
end
foundationFile = cfg.files.lowSpeedTemplate;
sourceFile = R5_StageMatForMatlab(sourceFile,'r5_template_bank');
S = load(sourceFile,'LowSpeedTemplateBank');
assert(isfield(S,'LowSpeedTemplateBank') && isfield(S.LowSpeedTemplateBank,'entry'), ...
    'R5:LowSpeedBankMissing','Low-speed template bank has no entry array.');
B = S.LowSpeedTemplateBank;
keep = false(numel(B.entry),1);
for i = 1:numel(B.entry)
    keep(i) = ismember(double(B.entry(i).sensorId),double(cfg.case.analysisSensors));
end
src = B.entry(keep);
blank = struct('blade_id',NaN,'sensor_id',NaN,'x_grid',[],'v_grid',[], ...
    'baseline',0,'valid_grid_mask',[],'domain_effective_mask',[], ...
    'sourceTemplateFile','','voltageUnit','V');
Template = struct('schema','R5_LOWSPEED_TEMPLATE_BANK_ADAPTER_V1', ...
    'sourceFile',sourceFile,'SensorBlade',repmat(blank,numel(src),1));
for i = 1:numel(src)
    a = src(i).sensorTemplate;
    x = double(a.x_grid(:));
    v = (double(a.v_grid(:))-double(a.baseline))/1000; % mV -> V
    valid = isfinite(x) & isfinite(v);
    if isfield(a,'x_domain')
        valid = valid & x >= double(a.x_domain(1)) & x <= double(a.x_domain(2));
    end
    Template.SensorBlade(i) = struct( ...
        'blade_id',double(src(i).bladeId), ...
        'sensor_id',double(src(i).sensorId), ...
        'x_grid',x,'v_grid',v,'baseline',0, ...
        'valid_grid_mask',valid,'domain_effective_mask',valid, ...
        'sourceTemplateFile',src(i).sourceTemplateFile, ...
        'voltageUnit','V');
end
% The target B4 template used by the dynamic foundation is more recent than
% the legacy bank.  Keep the bank for non-target anchors, but replace B4
% target waveforms with the foundation template so the dynamic baseline is
% exactly the one used by the corresponding WindowResult.
if isfile(foundationFile)
    FT = load(R5_StageMatForMatlab(foundationFile,'r5_target_template'),'Template');
    for j=1:numel(FT.Template.SensorBlade)
        a=FT.Template.SensorBlade(j); bid=double(a.blade_id); sid=double(a.sensor_id);
        ii=find([Template.SensorBlade.blade_id]==bid & [Template.SensorBlade.sensor_id]==sid,1);
        if isempty(ii), continue; end
        x=double(a.x_grid(:)); v=double(a.v_grid(:)); base=double(a.baseline);
        valid=isfinite(x)&isfinite(v);
        Template.SensorBlade(ii).x_grid=x;
        Template.SensorBlade(ii).v_grid=v;
        Template.SensorBlade(ii).baseline=base;
        Template.SensorBlade(ii).valid_grid_mask=valid;
        Template.SensorBlade(ii).domain_effective_mask=valid;
        Template.SensorBlade(ii).sourceTemplateFile=foundationFile;
    end
end
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'Template','-v7.3');
fprintf('20241106 R5 low-speed template adapter saved: %s\n',outputFile);
end
