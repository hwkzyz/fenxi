function T = R5_Localize_R4_Surface_FromLowSpeed(familyFile, templateFile, outputDir, b2GapMm, registrationFile)
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
load(familyFile,'F'); load(templateFile,'Template');
load(registrationFile,'G');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
rows={}; details=struct([]); k=0;
for i=1:numel(Template.SensorBlade)
    a=Template.SensorBlade(i); bid=double(a.blade_id); sid=double(a.sensor_id);
    xa=double(a.x_grid(:)); ya=(double(a.v_grid(:))-double(a.baseline))*1000;
    m=isfinite(xa)&isfinite(ya);
    if isfield(a,'valid_grid_mask'), m=m&logical(a.valid_grid_mask(:)); end
    if isfield(a,'domain_effective_mask'), m=m&logical(a.domain_effective_mask(:)); end
    if nnz(m)<20, continue; end
    ir=find([G.registration.sensor_id]==sid,1);
    assert(~isempty(ir),'R5:RegistrationMissing','No B2 registration for sensor %d.',sid);
    xLib=G.registration(ir).x_scale*(xa-G.registration(ir).tau_mm);
    x=F.x_mm; y=interp1(xLib(m),ya(m),x,'pchip',NaN);
    use=F.support_mask & isfinite(y);
    if nnz(use)<20, continue; end
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
    best=[inf NaN NaN];
    if bid==2
        gg=b2GapMm;
    end
    Xg=[ones(1,numel(gg)); 1./gg; log(gg/F.gref_mm)];
    for it=1:numel(tg)
        q=interp1(z,qord.',tg(it),'linear').';
        B=reshape(q,numel(x),3);
        yp=G.registration(ir).voltage_gain*(B*Xg) + G.registration(ir).voltage_offset_mv;
        e=sqrt(mean((yp(use,:)-y(use)).^2,1,'omitnan'));
        [emin,ig]=min(e);
        if emin<best(1), best=[emin gg(ig) tg(it)]; end
    end
    k=k+1; isB2=(bid==2);
    rows{k,1}=table(bid,sid,best(2),best(3),best(1),isB2,nnz(use)/numel(use), ...
        'VariableNames',{'blade_id','sensor_id','gap_or_effective_mm','latent_coordinate', ...
        'anchor_rmse_mv','b2_fixed_gap','support_fraction'});
    details(k).blade_id=bid; details(k).sensor_id=sid; details(k).x_mm=x;
    details(k).anchor_mv=y; details(k).support=use; details(k).gap=best(2); details(k).latent=best(3);
    q=interp1(z,qord.',best(3),'linear').'; details(k).B=reshape(q,numel(x),3);
end
T=vertcat(rows{:});
writetable(T,fullfile(outputDir,'R5_R4_surface_localization.csv'));
save(fullfile(outputDir,'R5_R4_surface_localization.mat'),'T','details','familyFile','templateFile','b2GapMm','registrationFile','-v7.3');
fprintf('R5 surface localization saved: %s\n',outputDir);
end
