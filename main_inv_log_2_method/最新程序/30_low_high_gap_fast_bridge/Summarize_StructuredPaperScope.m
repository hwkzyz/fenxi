function Summary=Summarize_StructuredPaperScope()
%SUMMARIZE_STRUCTUREDPAPERSCOPE Consolidate the 10/15/20 dB paper-scope runs.
root=fileparts(mfilename('fullpath'));outRoot=fullfile(root,'output');snrList=[10 15 20];
modes=["single_sync","single_async","dual_sync","sync_async"];
rows=repmat(row0(),numel(snrList)*numel(modes),1);at=0;
for snr=snrList
    for mode=modes
        at=at+1;folder=sprintf('%s_structured_%ddb_8rev',mode,snr);
        file=fullfile(outRoot,folder,sprintf('%s_structured.csv',mode));T=readtable(file);
        switch mode
            case "single_sync",ferr=abs(T.f_est-T.f_true);
            case "single_async",ferr=T.frequency_error;
            otherwise,ferr=T.max_frequency_error;
        end
        rows(at)=struct('mode',mode,'snr_db',snr,'case_count',height(T),...
            'success_count',nnz(T.success),'success_rate',mean(T.success),...
            'max_frequency_error_hz',max(ferr),'max_gap_error_mm',max(abs(T.g_est-T.g_true)),...
            'mean_elapsed_s',mean(T.elapsed_s),'max_elapsed_s',max(T.elapsed_s),...
            'used_samples',min(T.used_samples));
    end
end
Summary=struct2table(rows);out=fullfile(outRoot,'structured_paper_scope_summary');
if ~exist(out,'dir'),mkdir(out);end
writetable(Summary,fullfile(out,'structured_paper_scope_summary.csv'));
save(fullfile(out,'structured_paper_scope_summary.mat'),'Summary');disp(Summary);
if any(Summary.success_rate<1),error('Paper-scope structured regression is not fully successful.');end
end

function r=row0()
r=struct('mode',"",'snr_db',NaN,'case_count',NaN,'success_count',NaN,...
    'success_rate',NaN,'max_frequency_error_hz',NaN,'max_gap_error_mm',NaN,...
    'mean_elapsed_s',NaN,'max_elapsed_s',NaN,'used_samples',NaN);
end
