function outputFile = R5_Export_Main10_CompatibleLibrary(familyFile,localizationFile,baseLibraryFile,outputFile,targetBlade)
% Export a sidecar Main10 input without changing Main10 or Main07 files.
% The frozen sensor correction fields are retained; only the common static
% gap surface is replaced by the R4 anchor-localized waveform surface.
if nargin<1 || isempty(familyFile), error('R5:MissingFamily','familyFile is required.'); end
if nargin<2 || isempty(localizationFile), error('R5:MissingLocalization','localizationFile is required.'); end
if nargin<3 || isempty(baseLibraryFile), error('R5:MissingBaseLibrary','baseLibraryFile is required.'); end
if nargin<4 || isempty(outputFile), outputFile=fullfile(fileparts(mfilename('fullpath')),'results','R5_Main10_CompatibleLibrary.mat'); end
if nargin<5 || isempty(targetBlade), targetBlade=1; end
load(familyFile,'F'); load(localizationFile,'T'); S=load(baseLibraryFile,'CorrectedGapLibrary');
assert(isfield(F,'excludedTargetBlade') && double(F.excludedTargetBlade)==double(targetBlade), ...
    'R5:LOBOContract','Family does not exclude the exported target blade B%d.',targetBlade);
baseSensors = S.CorrectedGapLibrary.sensor;
use=T.blade_id==targetBlade;
assert(any(use),'R5:NoTargetAnchor','No low-speed anchors for B%d.',targetBlade);
zSharedSensor=double(T.latent_coordinate(use));
zSensor=zSharedSensor;
if ismember('sensorwise_latent_coordinate',T.Properties.VariableNames)
    zSensor=double(T.sensorwise_latent_coordinate(use));
end
zSensor=zSensor(isfinite(zSensor));
assert(numel(zSensor)>=2,'R5:InsufficientSensors', ...
    'At least two finite sensor anchors are required for B%d.',targetBlade);
zRange=range(zSensor); familyRange=range(F.latent_sorted);
zRangeFraction=zRange/max(familyRange,eps);
assert(zRangeFraction<=0.10,'R5:LatentSensorDisagreement', ...
    'B%d sensor latent range is %.2f%% of the calibrated family (>10%%).', ...
    targetBlade,100*zRangeFraction);
z=median(zSharedSensor,'omitnan');
assert(z>=min(F.latent_sorted) && z<=max(F.latent_sorted), ...
    'R5:LatentExtrapolation','Localized latent coordinate is outside the family range.');
[~,ord]=ismember(F.blade_order_sorted,F.blade_id);
q=interp1(F.latent_sorted,F.q(:,ord).',z,'linear').';
B=reshape(q,numel(F.x_mm),3);
R=S.CorrectedGapLibrary.responseSurface;
baseR=R;
assert(isequal(size(R.coeff),size(B)),'R5:SurfaceShape','Coefficient dimensions differ.');
assert(isequaln(double(R.xGrid(:)),double(F.x_mm(:))), ...
    'R5:XGridMismatch','Main10 and R5 surface x-grids differ.');
assert(isequaln(double(R.g0Mm),double(F.gref_mm)), ...
    'R5:ReferenceGapMismatch','Main10 and R5 reference gaps differ.');
assert(all(isfinite(B(:))),'R5:NonfiniteCoefficient','Export coefficients contain nonfinite values.');
R.coeff=B;
% Preserve all original library grids, sampled waveforms, and sensor-facing
% metadata. Main10 uses coeff as the forward static-gap surface; replacing
% calibration-support arrays would silently alter the library contract.
R.method='R5 PC1-ordered piecewise-linear coefficient-manifold sidecar for Main10 test';
S.CorrectedGapLibrary.responseSurface=R;
% The R5 sidecar is allowed to replace the response surface only.  All
% sensor-facing calibration/nuisance fields remain byte-for-byte identical.
assert(isequaln(S.CorrectedGapLibrary.sensor,baseSensors), ...
    'R5:SensorContractChanged','R5 sidecar changed frozen sensor calibration.');
S.CorrectedGapLibrary.r5=struct('family_file',familyFile,'localization_file',localizationFile, ...
    'target_blade',targetBlade,'localized_latent_coordinate',z, ...
    'sensor_latent_coordinates',zSensor(:).','sensor_latent_range',zRange, ...
    'sensor_latent_range_fraction',zRangeFraction,'sensor_consistency_limit_fraction',0.10, ...
    'note','Experimental latent coordinate; not physical tilt or absolute local clearance.');
CorrectedGapLibrary=S.CorrectedGapLibrary; %#ok<NASGU>
preservedFields={'xGrid','g0Mm','waveforms','effectiveWindow','dFdgGrid','dFdxGrid'};
for k=1:numel(preservedFields)
    name=preservedFields{k};
    assert(isequaln(CorrectedGapLibrary.responseSurface.(name),baseR.(name)), ...
        'R5:CalibrationArrayChanged','responseSurface.%s was altered.',name);
end
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'CorrectedGapLibrary','-v7.3');
fprintf('R5-compatible Main10 library: %s\n',outputFile);
end
