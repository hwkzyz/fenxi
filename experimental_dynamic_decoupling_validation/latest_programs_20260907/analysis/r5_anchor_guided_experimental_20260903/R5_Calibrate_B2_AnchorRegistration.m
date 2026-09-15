function G = R5_Calibrate_B2_AnchorRegistration(familyFile,templateFile,outputFile,b2GapMm,anchorBlade,sensorIds)
% Establish the OPR-x to static-library-x registration from B2 only.
% B2 is the sole blade with a known nominal gap.  This is a coordinate
% registration, not a clearance/tilt fit, and is fixed before R4 surface
% localization is applied to any other blade.
if nargin<1 || isempty(familyFile), familyFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_experimental_surface_family.mat'); end
if nargin<2 || isempty(templateFile), error('R5:MissingTemplate','A Main07 template is required.'); end
if nargin<3 || isempty(outputFile), outputFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_b2_anchor_registration.mat'); end
if nargin<4 || isempty(b2GapMm), b2GapMm=1.0; end
if nargin<5 || isempty(anchorBlade), anchorBlade=2; end
familyFileSource=familyFile; templateFileSource=templateFile;
familyFile=R5_StageMatForMatlab(familyFile,'r5_family');
templateFile=R5_StageMatForMatlab(templateFile,'r5_template');
load(familyFile,'F'); load(templateFile,'Template');
if nargin<6 || isempty(sensorIds), sensorIds=unique([Template.SensorBlade.sensor_id]); end
ib=find(F.blade_id==anchorBlade,1);
assert(~isempty(ib),'R5:AnchorMissing','Static family does not contain anchor B%d.',anchorBlade);
[~,ig]=min(abs(F.gap_mm-b2GapMm));
assert(abs(F.gap_mm(ig)-b2GapMm)<1e-9, ...
    'R5:B2GapUnavailable','The static family has no exact B2 calibration gap %.6g mm.',b2GapMm);
yStatic=F.waveforms_mv(:,ib,ig); xStatic=F.x_mm(:);
sensorId=double(sensorIds(:).');
row=cell(numel(sensorId),1);
registration=repmat(struct('sensor_id',NaN,'tau_mm',NaN,'x_scale',NaN, ...
    'voltage_gain',NaN,'voltage_offset_mv',NaN,'rmse_mv',NaN, ...
    'point_count',NaN,'b2_gap_mm',NaN),numel(sensorId),1);
for is=1:numel(sensorId)
    a=Template.SensorBlade([Template.SensorBlade.blade_id]==anchorBlade & [Template.SensorBlade.sensor_id]==sensorId(is));
    assert(numel(a)==1,'R5:AnchorMissing','B%d/S%d anchor missing.',anchorBlade,sensorId(is));
    x=double(a.x_grid(:)); y=(double(a.v_grid(:))-double(a.baseline))*1000;
    use=isfinite(x)&isfinite(y);
    if isfield(a,'valid_grid_mask'), use=use&logical(a.valid_grid_mask(:)); end
    if isfield(a,'domain_effective_mask'), use=use&logical(a.domain_effective_mask(:)); end
    obj=@(p) local_rmse(p,x(use),y(use),xStatic,yStatic);
    p=fminsearch(obj,[0,1],optimset('Display','off','TolX',1e-7,'TolFun',1e-5));
    p(1)=min(max(p(1),-1),1); p(2)=min(max(p(2),0.8),1.2);
    xLib=p(2)*(x-p(1)); yRaw=interp1(xStatic,yStatic,xLib,'pchip',NaN);
    fitUse=use&isfinite(yRaw);
    ab=[yRaw(fitUse),ones(nnz(fitUse),1)]\y(fitUse);
    yFit=ab(1)*yRaw+ab(2); rmse=sqrt(mean((yFit(fitUse)-y(fitUse)).^2));
    registration(is)=struct('sensor_id',sensorId(is),'tau_mm',p(1),'x_scale',p(2), ...
        'voltage_gain',ab(1),'voltage_offset_mv',ab(2),'rmse_mv',rmse, ...
        'point_count',nnz(fitUse),'b2_gap_mm',b2GapMm);
    row{is}=table(sensorId(is),p(1),p(2),ab(1),ab(2),rmse,nnz(fitUse), ...
        'VariableNames',{'sensor_id','tau_mm','x_scale','voltage_gain', ...
        'voltage_offset_mv','b2_anchor_rmse_mv','point_count'});
end
G=struct('schema','R5_ANCHOR_COORDINATE_REGISTRATION_V1', ...
    'family_file',familyFileSource,'template_file',templateFileSource,'b2_gap_mm',b2GapMm, ...
    'anchor_blade',anchorBlade,'registration',registration,'summary',vertcat(row{:}), ...
    'note','Known-gap anchor fixes spatial and affine-voltage registration before target localization.');
if ~exist(fileparts(outputFile),'dir'), mkdir(fileparts(outputFile)); end
save(outputFile,'G','-v7.3');
writetable(G.summary,[erase(outputFile,'.mat') '.csv']);
end

function e=local_rmse(p,x,y,xStatic,yStatic)
tau=min(max(p(1),-1),1); k=min(max(p(2),0.8),1.2);
yp=interp1(xStatic,yStatic,k*(x-tau),'pchip',NaN); use=isfinite(yp);
if nnz(use)<20, e=1e6; return; end
ab=[yp(use),ones(nnz(use),1)]\y(use); e=sqrt(mean((ab(1)*yp(use)+ab(2)-y(use)).^2));
end
