function file = R5_Plot_StaticWaveformComparison(x, anchor, truth, predicted, titleText, file)
%R5_PLOT_STATICWAVEFORMCOMPARISON Required method-level waveform evidence.
fig=figure('Visible','off','Color','w'); hold on; grid on; box on;
plot(x,anchor,'k--','LineWidth',1,'DisplayName','anchor');
plot(x,truth,'b-','LineWidth',1.4,'DisplayName','held-out measured');
plot(x,predicted,'r-','LineWidth',1.2,'DisplayName','R5 migration');
xlabel('x_{B2,OPR} (mm)'); ylabel('Voltage (mV)'); title(titleText); legend('Location','best');
exportgraphics(fig,file,'Resolution',180); close(fig);
end
