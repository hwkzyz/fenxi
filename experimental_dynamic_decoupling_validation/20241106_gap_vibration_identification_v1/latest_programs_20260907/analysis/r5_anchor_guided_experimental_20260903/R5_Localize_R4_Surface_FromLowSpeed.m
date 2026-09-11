function T = R5_Localize_R4_Surface_FromLowSpeed(familyFile, templateFile, outputDir, b2GapMm, registrationFile, anchorBlade, sensorIds, latentMode)
% Localize each blade response surface from OPR-aligned Main07 templates.
% B2 uses its known nominal gap; other blades use a joint (gap,latent) scan.
if nargin<1 || isempty(familyFile), familyFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_experimental_surface_family.mat'); end
if nargin<2 || isempty(templateFile)
    templateFile=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
        'inputs','prepared','foundation','step04_low_speed_template', ...
        'Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S123_20251222.mat');
end
if nargin<3 || isempty(outputDir), outputDir=fullfile(fileparts(mfilename('fullpath')),'results','surface_localization'); end
if nargin<4 || isempty(b2GapMm), b2GapMm=1.0; end
if nargin<5 || isempty(registrationFile), registrationFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_b2_anchor_registration.mat'); end
if nargin<6 || isempty(anchorBlade), anchorBlade=2; end
familyFileSource=familyFile; templateFileSource=templateFile; registrationFileSource=registrationFile;
familyFile=R5_StageMatForMatlab(familyFile,'r5_family'); templateFile=R5_StageMatForMatlab(templateFile,'r5_template'); registrationFile=R5_StageMatForMatlab(registrationFile,'r5_registration');
load(familyFile,'F'); load(templateFile,'Template');
load(registrationFile,'G');
if nargin<7 || isempty(sensorIds), sensorIds=unique([Template.SensorBlade.sensor_id]); end
if nargin<8 || isempty(latentMode), latentMode='sensor_conditioned'; end
if isscalar(b2GapMm), b2GapMm=repmat(double(b2GapMm),1,numel(sensorIds)); else, b2GapMm=double(b2GapMm(:).'); end
assert(numel(b2GapMm)==numel(sensorIds),'R5:B2GapSensorMismatch', ...
    'One anchor gap is required for each localization sensor.');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
rows={}; details=struct([]); k=0;
for i=1:numel(Template.SensorBlade)
    a=Template.SensorBlade(i); bid=double(a.blade_id); sid=double(a.sensor_id);
    if ~ismember(sid,double(sensorIds)), continue; end
    xa=double(a.x_grid(:)); ya=(double(a.v_grid(:))-double(a.baseline))*1000;
    m=isfinite(xa)&isfinite(ya);
    if isfield(a,'valid_grid_mask'), m=m&logical(a.valid_grid_mask(:)); end
    if isfield(a,'domain_effective_mask'), m=m&logical(a.domain_effective_mask(:)); end
    if nnz(m)<20, continue; end
    ir=find([G.registration.sensor_id]==sid,1);
    assert(~isempty(ir),'R5:RegistrationMissing','No B2 registration for sensor %d.',sid);
    x=F.x_mm;
    % The coefficient columns must follow the same order as the sorted
    % latent coordinate.  The previous ismember expression retained the
    % original blade order and paired every latent value with a different
    % surface, which invalidated the anchor localization.
    z=F.latent_sorted;
    [~,ord]=ismember(F.blade_order_sorted,F.blade_id);
    qord=F.q(:,ord);
    % R4 localizes on the observed waveform family.  Extrapolation beyond
    % that family creates artificial boundary minima and is not supported by
    % an experimental anchor waveform.
    tg=linspace(z(1),z(end),301);
    gg=linspace(min(F.gap_mm),max(F.gap_mm),281);
    best=[inf NaN NaN 0];
    if bid==anchorBlade
        gg=b2GapMm(sensorIds==sid);
    end
    offsets=0;
    if bid~=anchorBlade, offsets=linspace(-0.30,0.30,61); end
    Xg=[ones(1,numel(gg)); 1./gg; log(gg/F.gref_mm)];
    for io=1:numel(offsets)
        xLib=G.registration(ir).x_scale*(xa(m)-G.registration(ir).tau_mm-offsets(io));
        y=interp1(xLib,ya(m),x,'pchip',NaN);
        use=F.support_mask & isfinite(y);
        if nnz(use)<20, continue; end
        for it=1:numel(tg)
            q=interp1(z,qord.',tg(it),'linear').';
            B=reshape(q,numel(x),3);
            yp=G.registration(ir).voltage_gain*(B*Xg) + G.registration(ir).voltage_offset_mv;
            e=sqrt(mean((yp(use,:)-y(use)).^2,1,'omitnan'));
            [emin,ig]=min(e);
            if emin<best(1), best=[emin gg(ig) tg(it) offsets(io)]; end
        end
    end
    if ~isfinite(best(1)), continue; end
    % Store the latent/gap profile at the selected fixed target offset.
    xLib=G.registration(ir).x_scale*(xa(m)-G.registration(ir).tau_mm-best(4));
    y=interp1(xLib,ya(m),x,'pchip',NaN); use=F.support_mask & isfinite(y);
    profileRmse=nan(numel(tg),numel(gg));
    for it=1:numel(tg)
        q=interp1(z,qord.',tg(it),'linear').'; B=reshape(q,numel(x),3);
        yp=G.registration(ir).voltage_gain*(B*Xg)+G.registration(ir).voltage_offset_mv;
        profileRmse(it,:)=sqrt(mean((yp(use,:)-y(use)).^2,1,'omitnan'));
    end
    k=k+1; isB2=(bid==anchorBlade);
    rows{k,1}=table(bid,sid,best(2),best(3),best(1),isB2,nnz(use)/numel(use),best(4), ...
        'VariableNames',{'blade_id','sensor_id','gap_or_effective_mm','latent_coordinate', ...
        'anchor_rmse_mv','b2_fixed_gap','support_fraction','target_x_offset_mm'});
    details(k).blade_id=bid; details(k).sensor_id=sid; details(k).x_mm=x;
    details(k).anchor_mv=y; details(k).support=use; details(k).gap=best(2); details(k).latent=best(3); details(k).target_x_offset_mm=best(4);
    details(k).profile_z=tg; details(k).profile_rmse=profileRmse; details(k).profile_gap=gg;
    q=interp1(z,qord.',best(3),'linear').'; details(k).B=reshape(q,numel(x),3);
end
T=vertcat(rows{:});
T.sensorwise_latent_coordinate = T.latent_coordinate;
T.sensorwise_gap_or_effective_mm = T.gap_or_effective_mm;
T.latent_mode = repmat(string(latentMode), height(T), 1);
% Default: keep the response coordinate sensor-conditioned. A shared
% coordinate is an optional diagnostic/model reduction, not an automatic
% physical constraint imposed merely because blade_id is the same.
if strcmpi(latentMode, 'shared')
    for bid=unique(T.blade_id(:)).'
        idx=find(T.blade_id==bid);
        if numel(idx)<2, continue; end
        z0=details(idx(1)).profile_z(:); J=zeros(size(z0));
        for ii=1:numel(idx)
            p=details(idx(ii)).profile_rmse; [~,ig]=min(p,[],2);
            J=J+p(sub2ind(size(p),(1:numel(z0)).',ig)).^2;
        end
        [~,iz]=min(J); zShared=z0(iz);
        for ii=1:numel(idx)
            p=details(idx(ii)).profile_rmse; [rmse,ig]=min(p(iz,:));
            T.latent_coordinate(idx(ii))=zShared;
            T.gap_or_effective_mm(idx(ii))=details(idx(ii)).profile_gap(ig);
            T.anchor_rmse_mv(idx(ii))=rmse;
        end
    end
end
% Short filenames avoid the Windows MAX_PATH/HDF5 failure in the packaged
% experiment directory; the manifest remains the provenance index.
writetable(T,fullfile(outputDir,'loc.csv'));
familyFile=familyFileSource; templateFile=templateFileSource; registrationFile=registrationFileSource;
save(fullfile(outputDir,'loc.mat'),'T','details','familyFile','templateFile','b2GapMm','registrationFile','latentMode','-v7.3');
fprintf('R5 surface localization saved: %s\n',outputDir);
end
