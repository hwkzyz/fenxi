%% Overlay all calibration-library gap curves with low-speed rotor templates
clear; clc; close all;
root = fileparts(fileparts(mfilename('fullpath')));
libFile = fullfile(root,'results','Step06I_OffsetTiltShared_GapLibrary_20250527_B1_S136_nob_harddomain_fixed12pctroi_formal_20260830.mat');
S = load(libFile,'CorrectedGapLibrary'); C = S.CorrectedGapLibrary;
R = C.responseSurface; T = load(C.templateFile,'Template'); T = T.Template;
sids = [1 3 6];
SB = T.SensorBlade;
outDir = fullfile(root,'analysis','outputs','figures_all_gap_library_vs_low_speed'); if ~isfolder(outDir), mkdir(outDir); end
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 15]); tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
for ii=1:numel(sids)
    sid=sids(ii); cs=C.sensor(find(arrayfun(@(z) z.sensorId==sid,C.sensor),1)); ts=SB(find(arrayfun(@(z) z.sensor_id==sid,SB),1));
    xLow=ts.x_grid(:); v=ts.v_grid_baseline_removed(:)*1000;
    x=cs.x(:); xlib=cs.xLib(:); % library-supported low-speed ROI
    xq=R.xGrid(:); ug=R.trueGapMm(:); hold on;
    nexttile; hold on;
    cmap=parula(numel(ug));
    for j=1:numel(ug)
        y=cs.voltageGain*interp1(xq,R.YmeanAll(:,j),xlib,'linear',NaN);
        plot(x,y,'-','Color',0.65*cmap(j,:)+0.35,'LineWidth',0.45,'HandleVisibility','off');
    end
    plot(xLow,v,'k-','LineWidth',1.35,'DisplayName','rotor low-speed template');
    plot(cs.x,cs.vFitMv,'r--','LineWidth',1.15,'DisplayName',sprintf('fitted path (RMSE %.1f mV)',cs.lowFitRmseMv));
    [~,iMin] = min(cs.gEff); [~,iMax] = max(cs.gEff);
    plot(cs.x(iMin),cs.vFitMv(iMin),'ko','MarkerFaceColor',[0.1 0.1 0.1],...
        'DisplayName',sprintf('path min g %.3f mm',cs.gEff(iMin)));
    plot(cs.x(iMax),cs.vFitMv(iMax),'ks','MarkerFaceColor',[0.8 0.8 0.8],...
        'DisplayName',sprintf('path max g %.3f mm',cs.gEff(iMax)));
    xline(min(xq),'k:','library x-domain','HandleVisibility','off'); xline(max(xq),'k:','HandleVisibility','off');
    title(sprintf('20250527 CH%d | all %d gap curves, g = %.3f–%.3f mm',sid,numel(ug),min(ug),max(ug)),'Interpreter','none');
    xlabel('x (mm)'); ylabel('Voltage (mV)'); box on; set(gca,'TickDir','in','FontName','Times New Roman','FontSize',8.5); grid on;
    if ii==1
        legend('Location','northwest','Box','off');
        colormap(parula(numel(ug))); caxis([min(ug) max(ug)]); cb=colorbar; cb.Label.String='true gap g (mm)';
    end
end
exportgraphics(fig,fullfile(outDir,'AllGapLibrary_vs_LowSpeed_20250527_B1_S136.png'),'Resolution',300);
exportgraphics(fig,fullfile(outDir,'AllGapLibrary_vs_LowSpeed_20250527_B1_S136.pdf'),'ContentType','vector');
fprintf('Saved to %s\n',outDir);
