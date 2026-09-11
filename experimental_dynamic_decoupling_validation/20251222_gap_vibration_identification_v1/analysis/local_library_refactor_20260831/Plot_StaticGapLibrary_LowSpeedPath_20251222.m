% Visual audit of the static gap library and rotor low-speed paths.
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(scriptDir));
calFile = fullfile(rootDir,'results','prepared','calibration', ...
    'Step06I_OffsetTiltShared_GapLibrary_20251222_B1_S123_nob_formal_20260830.mat');
localFile = fullfile(scriptDir,'LocalNoVibrationGapLibrary_20251222_B1_S123_nob_20260831.mat');
outDir = fullfile(scriptDir,'figures_static_gap_library_20260901');
if ~exist(outDir,'dir'), mkdir(outDir); end
C = load(calFile).CorrectedGapLibrary.responseSurface;
L = load(localFile).LocalNoVibrationGapLibrary;
style = paper_style_local();
x = C.xGrid(:); gNom = C.trueGapMm(:); G = C.gEffByX;
nGap = numel(gNom); nBlade = numel(C.bladeIds);

% Figure 1: all mapped static paths, with the reference blade highlighted.
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 8.5]); hold on; box on;
cols=parula(nGap);
for ib=1:nBlade
    for j=1:nGap
        col=(ib-1)*nGap+j; if col>size(G,2), continue; end
        if ib==C.referenceBladeIndex, lw=1.3; ls='-'; else, lw=0.6; ls='--'; end
        plot(x,G(:,col),'Color',cols(j,:),'LineStyle',ls,'LineWidth',lw,'HandleVisibility','off');
    end
end
plot(nan,nan,'k-','LineWidth',1.3,'DisplayName',sprintf('Reference blade B%d',C.referenceBladeId));
plot(nan,nan,'k--','LineWidth',0.8,'DisplayName','Other blades');
colormap(cols); caxis([gNom(1) gNom(end)]); cb=colorbar; cb.Label.String='Nominal gap (mm)';
xlabel('x (mm)'); ylabel('Mapped gap path (mm)'); title('Static library: mapped gap paths');
legend('Location','eastoutside','Box','off'); format_axes_local(gca,style);
exportgraphics(fig,fullfile(outDir,'Fig1_StaticLibrary_AllBladeGapPaths.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig1_StaticLibrary_AllBladeGapPaths.pdf'),'ContentType','vector'); close(fig);

% Figure 2: measured low-speed template and localized reference path.
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 18]);
tiledlayout(numel(L.sensor),2,'TileSpacing','compact','Padding','compact');
for is=1:numel(L.sensor)
    s=L.sensor(is); xx=s.xGridMm(:);
    nexttile; plot(xx,s.T0MeasuredMv(:),'Color',[.12 .12 .12],'LineWidth',.85); box on;
    xlabel('x (mm)'); ylabel('T_0 (mV)'); title(sprintf('CH%d: low-speed template',s.sensorId)); format_axes_local(gca,style);
    nexttile; hold on; box on; yline(.90,'--','Color',[.55 .55 .55]); yline(2.30,'--','Color',[.55 .55 .55]);
    plot(xx,s.gRefPathMm(:),'Color',[.00 .28 .70],'LineWidth',1.15);
    xlabel('x (mm)'); ylabel('g_{ref}(x) (mm)'); title(sprintf('CH%d: localized path',s.sensorId)); ylim([.85 2.35]); format_axes_local(gca,style);
end
exportgraphics(fig,fullfile(outDir,'Fig2_LowSpeedWaveform_and_ReferencePaths.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig2_LowSpeedWaveform_and_ReferencePaths.pdf'),'ContentType','vector'); close(fig);

% Figure 3: low-speed paths relative to all reference-blade library curves.
fig=figure('Color','w','Units','centimeters','Position',[2 2 17 8.5]); hold on; box on;
for j=1:nGap, plot(x,G(:,(C.referenceBladeIndex-1)*nGap+j),'Color',[.78 .78 .78],'LineWidth',.55,'HandleVisibility','off'); end
sensorColors=[.85 .15 .12; 0 .45 .28; .55 .20 .70];
for is=1:numel(L.sensor), s=L.sensor(is); plot(s.xGridMm,s.gRefPathMm,'Color',sensorColors(is,:),'LineWidth',1.7,'DisplayName',sprintf('CH%d low-speed path',s.sensorId)); end
yline(.90,':','Color',[.35 .35 .35],'HandleVisibility','off'); yline(2.30,':','Color',[.35 .35 .35],'HandleVisibility','off');
xlabel('x (mm)'); ylabel('Gap (mm)'); title('Low-speed paths vs static library'); legend('Location','eastoutside','Box','off'); format_axes_local(gca,style);
exportgraphics(fig,fullfile(outDir,'Fig3_LowSpeedPaths_vs_StaticLibrary.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'Fig3_LowSpeedPaths_vs_StaticLibrary.pdf'),'ContentType','vector'); close(fig);

summary=table((1:numel(L.sensor)).',arrayfun(@(s)s.gCenterMm,L.sensor).',arrayfun(@(s)min(s.gRefPathMm),L.sensor).',arrayfun(@(s)max(s.gRefPathMm),L.sensor).','VariableNames',{'sensorId','gCenterMm','gRefMinMm','gRefMaxMm'});
writetable(summary,fullfile(outDir,'LowSpeedPath_Summary_20260901.csv')); disp(summary);

function style=paper_style_local(), style.fontName='Times New Roman'; style.tickFontSize=8.2; end
function format_axes_local(ax,style), set(ax,'FontName',style.fontName,'FontSize',style.tickFontSize,'Box','on','TickDir','in','LineWidth',.7,'Layer','top'); end
