function report = Plot_C2_FrequencyStabilityReferenceAudit(comparisonCsv,referenceCsv,outputDir,label)
%PLOT_C2_FREQUENCYSTABILITYREFERENCEAUDIT Audit primary-frequency stability.
% Foundation, V1 and R5 are plotted per window; the independent strain/BTT
% frequency is overlaid only as a posterior reference (never as an input).
if nargin<3||isempty(outputDir), outputDir=fullfile(fileparts(comparisonCsv),'frequency_audit'); end
if nargin<4||isempty(label), label='case'; end
if ~exist(outputDir,'dir'),mkdir(outputDir);end
T=readtable(comparisonCsv); n=height(T); w=T.window_id;
ref=nan(n,1);
if nargin>=2 && ~isempty(referenceCsv) && isfile(referenceCsv)
 R=readtable(referenceCsv); [~,ia,ib]=intersect(w,R.WindowID); ref(ia)=R.StrainFreq_Hz(ib);
end
names={'foundation_frequency_hz','v1_frequency_hz','r5_frequency_hz'}; labs={'Foundation','V1 gap-aware','R5'};
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 17 10]); hold on; box on;
co={[.45 .45 .45],[.1 .45 .75],[.85 .2 .15]};
for k=1:3, plot(w,T.(names{k}),'-o','Color',co{k},'MarkerSize',3,'LineWidth',1,'DisplayName',labs{k}); end
if any(isfinite(ref)), plot(w,ref,'k--','LineWidth',1.2,'DisplayName','independent strain/BTT'); end
xlabel('Window'); ylabel('Frequency (Hz)'); title([label ' | primary-frequency audit'],'FontWeight','normal'); legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',8,'TickDir','in','Box','on','XGrid','off','YGrid','off');
exportgraphics(fig,fullfile(outputDir,'C2_frequency_stability_reference.png'),'Resolution',180); close(fig);
M=table(w,T.foundation_frequency_hz,T.v1_frequency_hz,T.r5_frequency_hz,ref, ...
 'VariableNames',{'window_id','foundation_hz','v1_hz','r5_hz','reference_hz'});
M.r5_minus_reference_hz=M.r5_hz-M.reference_hz; M.v1_minus_reference_hz=M.v1_hz-M.reference_hz;
% Primary-frequency gate: this is an audit gate, not a frequency selector.
% It rejects out-of-band/branch-jump results but never overwrites a fitted value.
searchHz = [300 1000];
M.r5_in_search_band = isfinite(M.r5_hz) & M.r5_hz>=searchHz(1) & M.r5_hz<=searchHz(2);
M.r5_step_hz = [NaN; abs(diff(M.r5_hz))];
M.r5_continuity_pass = M.r5_in_search_band & (isnan(M.r5_step_hz) | M.r5_step_hz<=5);
% V1/Foundation and the independent reference are comparison evidence only.
% They must not veto R5: the methods intentionally model different physics,
% and the fixed-gap Foundation route may select a wrong branch by design.
M.r5_v1_agreement_pass = isfinite(M.r5_hz) & isfinite(M.v1_hz) & abs(M.r5_hz-M.v1_hz)<=2;
M.r5_reference_agreement_pass = ~isfinite(M.reference_hz) | (isfinite(M.r5_hz) & abs(M.r5_minus_reference_hz)<=2);
M.primary_frequency_gate = M.r5_continuity_pass;
writetable(M,fullfile(outputDir,'C2_frequency_stability_reference.csv'));
report=struct('label',label,'nWindows',n,'r5MedianHz',median(T.r5_frequency_hz,'omitnan'), ...
 'r5StdHz',std(T.r5_frequency_hz,'omitnan'),'referenceAvailable',any(isfinite(ref)), ...
 'searchHz',searchHz,'gatePass',all(M.primary_frequency_gate),'gatePassCount',sum(M.primary_frequency_gate), ...
 'outputDir',outputDir);
save(fullfile(outputDir,'C2_frequency_stability_reference.mat'),'report','M','-v7.3');
end
