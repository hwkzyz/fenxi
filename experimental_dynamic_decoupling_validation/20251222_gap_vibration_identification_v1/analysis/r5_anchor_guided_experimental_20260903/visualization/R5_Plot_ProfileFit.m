function file = R5_Plot_ProfileFit(selected, file)
%R5_PLOT_PROFILEFIT Plot measured and selected-state waveform in the same frame.
fig = figure('Visible','off','Color','w'); hold on; grid on; box on;
plot(selected.x_mm,selected.observed_mv,'k.','MarkerSize',5,'DisplayName','Main08 measured');
plot(selected.x_mm,selected.predicted_mv,'r-','LineWidth',1.2,'DisplayName','R5 selected state + vibration');
xlabel('x_{B2,OPR} (mm)'); ylabel('Voltage (mV)');
title(sprintf('R5 profile fit: %s=%.5g, EO=%.3g (%.1f Hz), A=%.4g mm', ...
    selected.stateAxisType,selected.stateValue,selected.EO,selected.frequencyHz,selected.amplitudeMm));
legend('Location','best'); exportgraphics(fig,file,'Resolution',180); close(fig);
end
