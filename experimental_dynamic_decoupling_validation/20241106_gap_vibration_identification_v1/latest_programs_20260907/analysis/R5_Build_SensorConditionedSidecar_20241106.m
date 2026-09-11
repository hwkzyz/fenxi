function L = R5_Build_SensorConditionedSidecar_20241106(familyFile, localizationFile, registrationFile, outputFile)
%R5_BUILD_SENSORCONDITIONEDSIDECAR_20241106 Store one response state per sensor.
% This is not a Main10-compatible library. It preserves z_{s,b} explicitly.
cfg = Config_20241106();
if nargin < 1 || isempty(familyFile), familyFile=fullfile(cfg.paths.results,'r5_experimental_surface_family.mat'); end
if nargin < 2 || isempty(localizationFile), localizationFile=fullfile(cfg.paths.results,'r5_surface_localization','loc.mat'); end
if nargin < 3 || isempty(registrationFile), registrationFile=fullfile(cfg.paths.results,'r5_b2_anchor_registration.mat'); end
if nargin < 4 || isempty(outputFile), outputFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
F=load(R5_StageMatForMatlab(familyFile,'r5_family'),'F'); X=load(R5_StageMatForMatlab(localizationFile,'r5_localization'),'T','details','latentMode'); G=load(R5_StageMatForMatlab(registrationFile,'r5_registration'),'G'); muBySensor=load_static_mu_local(cfg.files.gapLibrary,cfg.case.targetBlade);
assert(strcmpi(X.latentMode,'sensor_conditioned'),'R5:LatentMode','Sidecar requires sensor-conditioned localization.');
sidecar=struct('schema','R5_SENSOR_CONDITIONED_SIDECAR_V1','family_file',familyFile,...
    'localization_file',localizationFile,'registration_file',registrationFile,...
    'latent_mode',X.latentMode,'x_mm',F.F.x_mm,'gap_mm',F.F.gap_mm,'gref_mm',F.F.gref_mm,...
    'sensor',{{}});
for sid=unique(X.T.sensor_id(:)).'
    idx=find(X.T.sensor_id==sid); sensor=struct('sensor_id',sid,'rows',X.T(idx,:),... 
        'registration',G.G.registration([G.G.registration.sensor_id]==sid),...
        'coeff_by_state',{{}},...
        'state',repmat(struct('blade_id',NaN,'z',NaN,'gap_mm',NaN,'rmse_mv',NaN,'B',[],...
        'mu_gap_per_x_mm',0,'gap_tau_mm',NaN,'target_x_offset_mm',0),0,1));
    for j=1:numel(idx)
        d=X.details(idx(j));
        muState=0; tauState=NaN;
        if isfield(d,'mu_gap_per_x_mm'), muState=double(d.mu_gap_per_x_mm); end
        if isfield(d,'gap_tau_mm'), tauState=double(d.gap_tau_mm); end
        sensor.state(j)=struct('blade_id',X.T.blade_id(idx(j)),... 
            'z',X.T.sensorwise_latent_coordinate(idx(j)),...
            'gap_mm',X.T.sensorwise_gap_or_effective_mm(idx(j)),...
            'rmse_mv',X.T.anchor_rmse_mv(idx(j)),'B',d.B,...
            'mu_gap_per_x_mm',muState,'gap_tau_mm',tauState,'target_x_offset_mm',get_detail_offset_local(d));
        sensor.coeff_by_state{j}=d.B;
    end
    km=find([muBySensor.sensorId]==sid,1);
    if ~isempty(km)
        sensor.registration.mu_gap_per_x_mm=double(muBySensor(km).muGapPerXMm);
        if isfield(muBySensor,'tauMm'), sensor.registration.gap_tau_mm=double(muBySensor(km).tauMm); end
        for jj=1:numel(sensor.state)
            sensor.state(jj).mu_gap_per_x_mm=sensor.registration.mu_gap_per_x_mm;
            if isfield(sensor.registration,'gap_tau_mm'), sensor.state(jj).gap_tau_mm=sensor.registration.gap_tau_mm; end
        end
    else, sensor.registration.mu_gap_per_x_mm=0; end
    sidecar.sensor{end+1}=sensor; %#ok<AGROW>
end
L=struct('SensorConditionedLibrary',sidecar);
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'-struct','L','-v7.3');
end

function v=get_detail_offset_local(d)
v=0; if isfield(d,'target_x_offset_mm') && isfinite(d.target_x_offset_mm), v=double(d.target_x_offset_mm); end
end

function S=load_static_mu_local(file,targetBlade)
S=repmat(struct('sensorId',NaN,'muGapPerXMm',0,'tauMm',NaN),0,1);
if ~isfile(file), return; end
Q=load(R5_StageMatForMatlab(file,'r5_gap_mu'));
if isfield(Q,'CorrectedGapLibrary') && isfield(Q.CorrectedGapLibrary,'sensor')
    S=Q.CorrectedGapLibrary.sensor;
elseif isfield(Q,'GapCalibrationBank') && isfield(Q.GapCalibrationBank,'entry')
    E=Q.GapCalibrationBank.entry; keep=[E.bladeId]==targetBlade; E=E(keep); S=repmat(struct('sensorId',NaN,'muGapPerXMm',0,'tauMm',NaN),numel(E),1);
    for k=1:numel(E), S(k).sensorId=double(E(k).sensorId); S(k).muGapPerXMm=double(E(k).sensor.muGapPerXMm); if isfield(E(k).sensor,'tauMm'), S(k).tauMm=double(E(k).sensor.tauMm); end; end
end
end
