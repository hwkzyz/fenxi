function A = Audit_Blind_Fig02_Fig03_Data()
%AUDIT_BLIND_FIG02_FIG03_DATA Check numerical integrity before plotting.
root=fileparts(mfilename('fullpath')); d=fullfile(root,'paper_figures','updated_source_data','formal_blind_batch');
f2=readtable(fullfile(d,'R2_formal_blind.csv')); f3=readtable(fullfile(d,'R3_formal_blind.csv')); f4=readtable(fullfile(d,'R4_formal_blind.csv'));
assert(all(isfinite(f2.eo_true)) && all(isfinite(f2.eo_hat))); assert(all(f2.eo_hat>=1 & f2.eo_hat<=30));
assert(all(isfinite(f3.eo_true)) && all(isfinite(f3.eo_hat))); assert(all(f3.eo_hat>=1 & f3.eo_hat<=30));
assert(all(isfinite(str2double(string(f4.eo_true))))); assert(all(isfinite(str2double(string(f4.eo_hat)))));
S=table;
S.protocol=["R2";"R3";"R4"];
S.n=[height(f2);height(f3);height(f4)];
S.eo_recovery=[mean(f2.eo_true==f2.eo_hat);mean(f3.eo_true==f3.eo_hat);mean(string(f4.eo_true)==string(f4.eo_hat))];
S.identified=[mean(string(f2.status)=="identified");mean(string(f3.status)=="identified");mean(string(f4.identifiability_status)=="identified")];
S.max_abs_gap_error_mm=[max(abs(f2.gap_error_mm));max(abs(f3.gap_error_mm));max(abs(f4.g_high_error_mm))];
S.max_abs_amp_error_mm=[max(abs(f2.A_hat_mm-f2.A_true_mm));max(abs(f3.A_hat_mm-f3.A_true_mm));max(abs(str2double(string(f4.A_hat))-str2double(string(f4.A_true))) )];
S.rmse_median_V=[median(f2.rmse_V);median(f3.rmse_V);median(f4.high_rmse_V)];
out=fullfile(root,'paper_figures','updated_source_data','fig02_fig03_blind_review'); if ~exist(out,'dir'),mkdir(out);end
writetable(S,fullfile(out,'blind_integrity_summary.csv')); writetable(S,fullfile(out,'blind_integrity_summary.md'),'FileType','text');
fid=fopen(fullfile(out,'blind_integrity_summary.md'),'w'); fprintf(fid,'# Blind data integrity summary\n\n'); fprintf(fid,'| Protocol | n | EO recovery | identified | max gap error (mm) | max amplitude error (mm) | median RMSE (V) |\n|---|---:|---:|---:|---:|---:|---:|\n'); for i=1:3, fprintf(fid,'| %s | %d | %.3f | %.3f | %.6g | %.6g | %.6g |\n',S.protocol(i),S.n(i),S.eo_recovery(i),S.identified(i),S.max_abs_gap_error_mm(i),S.max_abs_amp_error_mm(i),S.rmse_median_V(i)); end; fclose(fid);
A=struct('summary',S,'outputDir',out); disp(S);
end
