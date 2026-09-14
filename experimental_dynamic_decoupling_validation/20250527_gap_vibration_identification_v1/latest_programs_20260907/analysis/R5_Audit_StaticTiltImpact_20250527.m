function report = R5_Audit_StaticTiltImpact_20250527(sidecarFile, templateFile, outputDir)
%R5_AUDIT_STATICTILTIMPACT_20250527 Quantify frozen static-slope influence.
% This diagnostic never changes the production R5 model or fitted results.
cfg = Config_20250527();
if nargin<1 || isempty(sidecarFile), sidecarFile=fullfile(cfg.paths.results,'r5_sensor_conditioned_sidecar.mat'); end
if nargin<2 || isempty(templateFile), templateFile=cfg.files.lowSpeedTemplate; end
if nargin<3 || isempty(outputDir), outputDir=fullfile(cfg.paths.results,'static_tilt_impact_audit'); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

S=load(R5_StageMatForMatlab(sidecarFile,'r5_tilt_sidecar'),'SensorConditionedLibrary');
Q=load(R5_StageMatForMatlab(templateFile,'r5_tilt_template'),'Template');
L=S.SensorConditionedLibrary; Template=Q.Template; dgProbe=0.02;
n=numel(L.sensor); rows=repmat(struct('sensor_id',NaN,'gap_reference_mm',NaN, ...
    'mu_gap_per_x_mm',NaN,'gap_min_mm',NaN,'gap_max_mm',NaN, ...
    'tilted_gap_min_mm',NaN,'tilted_gap_max_mm',NaN,'tilted_support_fraction',NaN, ...
    'increment_relative_rms_difference',NaN,'increment_correlation',NaN),n,1);

fig=figure('Visible','off','Color','w','Units','centimeters','Position',[1 1 17 15.8]);
tl=tiledlayout(fig,n,2,'TileSpacing','compact','Padding','compact');
for i=1:n
    sc=L.sensor{i}; sid=double(sc.sensor_id);
    js=find([sc.state.blade_id]==cfg.case.targetBlade,1);
    jt=find([Template.SensorBlade.sensor_id]==sid & ...
        [Template.SensorBlade.blade_id]==cfg.case.targetBlade,1);
    assert(~isempty(js)&&~isempty(jt),'R5:TiltAuditInput','Missing B%d/S%d state or template.',cfg.case.targetBlade,sid);
    st=sc.state(js); reg=sc.registration; tpl=Template.SensorBlade(jt);
    x=double(tpl.x_grid(:)); xOffset=0;
    if isfield(st,'target_x_offset_mm') && isfinite(st.target_x_offset_mm), xOffset=double(st.target_x_offset_mm); end
    xr=x-xOffset; xf=double(reg.x_scale).*(xr-double(reg.tau_mm));
    g0=double(st.gap_mm); mu=double(st.mu_gap_per_x_mm); tau=double(st.gap_tau_mm);
    gTilt=g0+mu.*(xr-tau); gMin=min(L.gap_mm); gMax=max(L.gap_mm);
    support=xf>=min(L.x_mm)&xf<=max(L.x_mm)&gTilt>=gMin&gTilt+dgProbe<=gMax;
    dScalar=eval_surface(st.B,L.x_mm,L.gref_mm,repmat(g0+dgProbe,size(xf)),xf)- ...
        eval_surface(st.B,L.x_mm,L.gref_mm,repmat(g0,size(xf)),xf);
    dTilt=nan(size(x));
    dTilt(support)=eval_surface(st.B,L.x_mm,L.gref_mm,gTilt(support)+dgProbe,xf(support))- ...
        eval_surface(st.B,L.x_mm,L.gref_mm,gTilt(support),xf(support));
    q=support&isfinite(dScalar)&isfinite(dTilt);
    rel=rms(dTilt(q)-dScalar(q))/max(rms(dScalar(q)),eps);
    cc=corrcoef(dTilt(q),dScalar(q)); if numel(cc)<4, rho=NaN; else, rho=cc(1,2); end
    rows(i)=struct('sensor_id',sid,'gap_reference_mm',g0,'mu_gap_per_x_mm',mu, ...
        'gap_min_mm',gMin,'gap_max_mm',gMax,'tilted_gap_min_mm',min(gTilt), ...
        'tilted_gap_max_mm',max(gTilt),'tilted_support_fraction',mean(support), ...
        'increment_relative_rms_difference',rel,'increment_correlation',rho);

    ax=nexttile(tl,2*i-1); hold(ax,'on');
    patch(ax,[min(x) max(x) max(x) min(x)],[gMin gMin gMax gMax],[.92 .92 .92], ...
        'EdgeColor','none','DisplayName','training range');
    plot(ax,x,g0+zeros(size(x)),'--','Color',[.20 .45 .72],'LineWidth',1.1,'DisplayName','scalar baseline');
    plot(ax,x,gTilt,'Color',[.82 .24 .18],'LineWidth',1.1,'DisplayName','frozen slope');
    ylabel(ax,'Gap (mm)'); title(ax,sprintf('S%d: static gap path',sid),'FontWeight','normal');
    if i==1, legend(ax,'Location','best','FontSize',8); end

    ax=nexttile(tl,2*i); hold(ax,'on');
    plot(ax,x,dScalar,'--','Color',[.20 .45 .72],'LineWidth',1.1,'DisplayName','scalar baseline');
    plot(ax,x,dTilt,'Color',[.82 .24 .18],'LineWidth',1.1,'DisplayName','frozen slope');
    ylabel(ax,'Voltage increment (mV)'); title(ax,sprintf('S%d: +%.2f mm response',sid,dgProbe),'FontWeight','normal');
    if i==1, legend(ax,'Location','best','FontSize',8); end
end
for k=1:2*n
    ax=nexttile(tl,k); xlabel(ax,'x (mm)'); box(ax,'on');
    set(ax,'FontName','Times New Roman','FontSize',8,'TickDir','in','LineWidth',.75,'XGrid','off','YGrid','off');
end
exportgraphics(fig,fullfile(outputDir,'A3_static_tilt_impact.png'),'Resolution',300);
exportgraphics(fig,fullfile(outputDir,'A3_static_tilt_impact.pdf'),'ContentType','vector'); close(fig);
T=struct2table(rows); writetable(T,fullfile(outputDir,'A3_static_tilt_impact.csv'));
report=struct('schema','R5_STATIC_TILT_IMPACT_AUDIT_V1','targetBlade',cfg.case.targetBlade, ...
    'dgProbeMm',dgProbe,'rows',T,'outputDir',outputDir);
save(fullfile(outputDir,'A3_StaticTiltImpactAudit.mat'),'report','T','-v7.3');
end

function y=eval_surface(B,xGrid,gref,g,xq)
y=zeros(size(xq));
for k=1:3
    bk=interp1(xGrid(:),B(:,k),xq,'linear',NaN);
    if k==1, y=y+bk; elseif k==2, y=y+bk./g; else, y=y+bk.*log(g/gref); end
end
end
