function T=R5_SensorGate_Sensitivity(localizationFile,outputFile,familyFile)
% Sensitivity of the declared normalized three-sensor latent gate.
if nargin<1||isempty(localizationFile), localizationFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_frontend_only','surface_localization','R5_R4_surface_localization.csv'); end
if nargin<2||isempty(outputFile), outputFile=fullfile(fileparts(mfilename('fullpath')),'results','r5_diagnostics_20260904','R5_SensorGate_Sensitivity.csv'); end
if nargin<3||isempty(familyFile), familyFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_frontend_only','r4_experimental_surface_family.mat'); end
A=readtable(localizationFile); L=load(familyFile,'F'); denom=max(range(L.F.latent_sorted),eps);
ids=unique(A.blade_id); spread=nan(numel(ids),1); nSensors=zeros(numel(ids),1);
for k=1:numel(ids), z=A.latent_coordinate(A.blade_id==ids(k)); z=z(isfinite(z)); nSensors(k)=numel(z); if ~isempty(z), spread(k)=range(z)/denom; end, end
threshold=[.05 .075 .10 .125 .15].'; passCount=zeros(size(threshold));
for k=1:numel(threshold), passCount(k)=sum(spread<=threshold(k)); end
T=table(threshold,passCount,repmat(numel(ids),numel(threshold),1), ...
    'VariableNames',{'threshold_fraction','passing_blades','total_blades'});
if ~isfolder(fileparts(outputFile)), mkdir(fileparts(outputFile)); end
writetable(T,outputFile); writetable(table(ids,spread,nSensors,'VariableNames',{'blade_id','normalized_spread','sensor_count'}), ...
    fullfile(fileparts(outputFile),'R5_SensorGate_ByBlade.csv'));
end
