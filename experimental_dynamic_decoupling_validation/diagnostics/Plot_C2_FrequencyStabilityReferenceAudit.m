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
writetable(M,fullfile(outputDir,'C2_frequency_stability_reference.csv'));
report=struct('label',label,'nWindows',n,'r5MedianHz',median(T.r5_frequency_hz,'omitnan'), ...
 'r5StdHz',std(T.r5_frequency_hz,'omitnan'),'referenceAvailable',any(isfinite(ref)), ...
 'outputDir',outputDir);
save(fullfile(outputDir,'C2_frequency_stability_reference.mat'),'report','M','-v7.3');
end
