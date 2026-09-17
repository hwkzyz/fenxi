function Result=Summarize_StructuredGapModelAblation(snrDb,seedOffset,truthMode)
%SUMMARIZE_STRUCTUREDGAPMODELABLATION Combine checkpointed method runs.
if nargin<1,snrDb=15;end
if nargin<2,seedOffset=0;end
if nargin<3,truthMode="absolute";end
truthMode=string(truthMode);
root=fileparts(mfilename('fullpath'));outputRoot=fullfile(root,'output');
methods=["fixed_low_gap","absolute_joint","low_template_increment"];
allRows=cell(numel(methods),1);
for i=1:numel(methods)
    folder=sprintf('structured_gap_ablation_%gdb_seed%d_%s_method_%s',...
        snrDb,seedOffset,truthMode,methods(i));
    file=fullfile(outputRoot,folder,'structured_gap_ablation.mat');
    if ~isfile(file),error('Missing ablation checkpoint: %s',file);end
    data=load(file,'Audit');allRows{i}=data.Audit;
end
Audit=vertcat(allRows{:});
Summary=repmat(struct('method',"",'case_count',0,'success_count',0,...
    'success_rate',NaN,'max_frequency_error_hz',NaN,'max_gap_error_mm',NaN,...
    'max_amplitude_relative_error',NaN,'median_voltage_rmse',NaN,...
    'median_fit_elapsed_s',NaN),numel(methods),1);
for i=1:numel(methods)
    q=Audit(Audit.method==methods(i),:);
    Summary(i)=struct('method',methods(i),'case_count',height(q),...
        'success_count',nnz(q.success),'success_rate',mean(q.success),...
        'max_frequency_error_hz',max(q.frequency_error_hz),...
        'max_gap_error_mm',max(q.gap_error_mm),...
        'max_amplitude_relative_error',max(q.max_amplitude_relative_error),...
        'median_voltage_rmse',median(q.rmse),...
        'median_fit_elapsed_s',median(q.elapsed_s));
end
Summary=struct2table(Summary);
out=fullfile(outputRoot,sprintf('structured_gap_ablation_summary_%gdb_seed%d_%s',...
    snrDb,seedOffset,truthMode));if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'structured_gap_ablation_combined.csv'));
writetable(Summary,fullfile(out,'structured_gap_ablation_method_summary.csv'));
save(fullfile(out,'structured_gap_ablation_summary.mat'),'Audit','Summary');
disp(Audit);disp(Summary);Result=struct('Audit',Audit,'Summary',Summary,'outputDir',out);
end
