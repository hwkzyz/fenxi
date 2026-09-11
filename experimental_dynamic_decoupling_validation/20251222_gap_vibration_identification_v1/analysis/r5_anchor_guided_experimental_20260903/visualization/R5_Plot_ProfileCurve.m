function file = R5_Plot_ProfileCurve(profile, file)
%R5_PLOT_PROFILECURVE Per-window J(state) evidence for joint selection.
P=profile.stateProfiles; good=isfinite([P.J]);
fig=figure('Visible','off','Color','w'); hold on; grid on; box on;
plot([P(good).stateValue],[P(good).J],'o-','LineWidth',1.2,'MarkerSize',4);
if isfinite(profile.selected.stateValue)
    xline(profile.selected.stateValue,'r--','LineWidth',1.1,'DisplayName','selected state');
end
xlabel('state axis value'); ylabel('profile SSE (mV^2)'); title('R5 profile objective by candidate state');
exportgraphics(fig,file,'Resolution',180); close(fig);
end
