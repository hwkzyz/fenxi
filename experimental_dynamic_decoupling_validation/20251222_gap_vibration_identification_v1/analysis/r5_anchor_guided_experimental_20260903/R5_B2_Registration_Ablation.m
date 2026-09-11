function T = R5_B2_Registration_Ablation(familyFile,templateFile,registrationFile,baseLibraryFile,outputDir)
% Apples-to-apples B2 registration comparison on the same low-speed anchors.
% R5 and frozen Main10 use the same B2 static family member, support and mV data;
% only their coordinate/voltage mappings differ.
if nargin<1||isempty(familyFile), familyFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_frontend_only','r4_experimental_surface_family.mat'); end
if nargin<2||isempty(templateFile), templateFile=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'inputs','prepared','foundation','step04_low_speed_template','Template_OPRCenterStd_AdaptiveSG_LowSpeed_AllBlades_S123_20251222.mat'); end
if nargin<3||isempty(registrationFile), registrationFile=fullfile(fileparts(mfilename('fullpath')),'results','r4_frontend_only','r4_b2_anchor_registration.mat'); end
if nargin<4||isempty(baseLibraryFile), baseLibraryFile=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'inputs','calibration','Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123.mat'); end
if nargin<5||isempty(outputDir), outputDir=fullfile(fileparts(mfilename('fullpath')),'results','registration_ablation'); end
if ~isfolder(outputDir), mkdir(outputDir); end
load(familyFile,'F'); load(templateFile,'Template'); load(registrationFile,'G'); load(baseLibraryFile,'CorrectedGapLibrary');
S=CorrectedGapLibrary.sensor; xFamily=F.x_mm(:); ig=find(abs(F.gap_mm-1)<1e-12,1);
assert(~isempty(ig),'R5:ReferenceGapMissing','Family has no 1.0-mm reference gap.');
ib=find(F.blade_id==2,1); yStatic=F.waveforms_mv(:,ib,ig);
rows=table();
for k=1:numel(G.registration)
    r=G.registration(k); a=Template.SensorBlade([Template.SensorBlade.blade_id]==2 & [Template.SensorBlade.sensor_id]==r.sensor_id);
    assert(~isempty(a),'R5:B2TemplateMissing','B2 template missing for CH%d.',r.sensor_id); a=a(1);
    xr=double(a.x_grid(:)); yr=(double(a.v_grid(:))-double(a.baseline))*1000;
    m=isfinite(xr)&isfinite(yr); if isfield(a,'valid_grid_mask'), m=m&logical(a.valid_grid_mask(:)); end; if isfield(a,'domain_effective_mask'), m=m&logical(a.domain_effective_mask(:)); end
    % Forward-map the same raw sensor points through each registration into
    % the common static x coordinate, then compare in the measured mV domain.
    j=find([S.sensorId]==r.sensor_id,1); s=S(j); xRaw=xr(m); yy=yr(m);
    xRfamily=r.x_scale*(xRaw-r.tau_mm); xMfamily=s.xScale*(xRaw-s.tauMm);
    rawR=interp1(xFamily,yStatic,xRfamily,'pchip',NaN); rawM=interp1(xFamily,yStatic,xMfamily,'pchip',NaN);
    predR=r.voltage_gain*rawR+r.voltage_offset_mv; predM=s.voltageGain*rawM+s.voltageOffsetMv;
    valid=isfinite(yy)&isfinite(predR)&isfinite(predM);
    yy=yy(valid); predR=predR(valid); predM=predM(valid); x=xRaw(valid);
    % Main10's actual dg=0 branch returns vLow directly.  Its baseline
    % residual is therefore the template self-residual, not a second absolute
    % registration fit.  The direct frozen-map curve is retained only as a
    % diagnostic and is explicitly not a runtime-comparable A/B result.
    main10BaselineRmse=sqrt(mean((yy-yy).^2));
    rows=[rows; table(r.sensor_id,nnz(valid),sqrt(mean((predR-yy).^2)),main10BaselineRmse, ...
        sqrt(mean((predR-predM).^2)),sqrt(mean((rawR(valid)-rawM(valid)).^2)), ...
        'VariableNames',{'sensor_id','point_count','r5_anchor_rmse_mv','main10_runtime_baseline_rmse_mv','r5_vs_direct_frozen_map_rmse_mv','raw_surface_map_rmse_mv'})]; %#ok<AGROW>
    f=figure('Visible','off','Color','w'); plot(x,yy,'k','DisplayName','same B2 low-speed anchor'); hold on; grid on; box on; plot(x,predR,'r','DisplayName','R5 registration'); plot(x,yy,'b--','DisplayName','Main10 runtime dg=0 (vLow)'); plot(x,predM,'Color',[0.3 0.3 0.3],'LineStyle',':','DisplayName','direct frozen-map diagnostic only'); xlabel('x / mm'); ylabel('mV'); title(sprintf('B2 registration contract audit, CH%d',r.sensor_id)); legend('Location','best'); exportgraphics(f,fullfile(outputDir,sprintf('R5_B2_Registration_AB_CH%d.png',r.sensor_id)),'Resolution',180); close(f);
end
writetable(rows,fullfile(outputDir,'R5_B2_Registration_Ablation.csv')); T=rows;
end
