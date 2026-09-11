function file = R5_Plot_CandidateStates(x, anchorY, candidates, selected, file)
%R5_PLOT_CANDIDATESTATES Save an auditable anchor/candidate waveform comparison.
if nargin < 5 || isempty(file), file = fullfile(pwd,'R5_candidate_waveforms.png'); end
fig = figure('Visible','off','Color','w'); hold on; grid on; box on;
plot(x, anchorY, 'k-', 'LineWidth', 1.5, 'DisplayName','anchor');
for k = 1:numel(candidates)
    if isempty(candidates(k).waveform_mv), continue; end
    good = strcmp(candidates(k).status,'candidate');
    if isstruct(selected) && isfinite(selected.stateValue) && candidates(k).stateValue == selected.stateValue
        sty = 'r-'; lw = 1.8;
    elseif good
        sty = '-'; lw = 0.8;
    else
        sty = ':'; lw = 0.5;
    end
    plot(x, candidates(k).waveform_mv, sty, 'LineWidth',lw, ...
        'DisplayName',sprintf('%s=%.5g',candidates(k).stateAxisType,candidates(k).stateValue));
end
xlabel('x_{B2,OPR} (mm)'); ylabel('Voltage (mV)'); title('R5 anchor-guided candidate waveforms');
legend('Location','bestoutside'); exportgraphics(fig,file,'Resolution',180); close(fig);
end
