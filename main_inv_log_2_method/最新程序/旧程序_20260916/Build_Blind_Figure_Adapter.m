function T = Build_Blind_Figure_Adapter()
%BUILD_BLIND_FIGURE_ADAPTER Create schema-stable panel summaries.
% Missing phase/noise-grid values remain NaN; no quantities are synthesized.
root=fileparts(mfilename('fullpath'));
src=fullfile(root,'paper_figures','updated_source_data','formal_blind_batch');
out=fullfile(root,'paper_figures','updated_source_data','fig02_fig03_blind_review');
if ~exist(out,'dir'), mkdir(out); end
r2=readtable(fullfile(src,'R2_formal_blind.csv')); r3=readtable(fullfile(src,'R3_formal_blind.csv')); r4=readtable(fullfile(src,'R4_formal_blind.csv'));
n=height(r2)+height(r3)+height(r4);
T=table('Size',[n 15],'VariableTypes',repmat({'string'},1,15),'VariableNames',{'figure_panel','protocol','case_id','eo_true','eo_hat','frequency_true_hz','frequency_hat_hz','gap_true_mm','gap_hat_mm','amplitude_true_mm','amplitude_hat_mm','phase_true_rad','phase_hat_rad','rmse_V','status'});
i=0;
for k=1:height(r2), i=i+1; T(i,:)=mk("Fig2d-i","R2",r2.case_id(k),r2.eo_true(k),r2.eo_hat(k),r2.gap_true_mm(k),r2.gap_hat_mm(k),r2.A_true_mm(k),r2.A_hat_mm(k),r2.rmse_V(k),r2.status(k)); end
for k=1:height(r3), i=i+1; T(i,:)=mk("Fig2d-i","R3",r3.case_id(k),r3.eo_true(k),r3.eo_hat(k),r3.gap_true_mm(k),r3.gap_hat_mm(k),r3.A_true_mm(k),r3.A_hat_mm(k),r3.rmse_V(k),r3.status(k)); end
for k=1:height(r4), i=i+1; T(i,:)=mk("Fig3f","R4",r4.case_id(k),str2double(string(r4.eo_true(k))),str2double(string(r4.eo_hat(k))),r4.g_high_true_mm(k),r4.g_high_hat_mm(k),str2double(string(r4.A_true(k))),str2double(string(r4.A_hat(k))),r4.high_rmse_V(k),r4.identifiability_status(k)); end
writetable(T,fullfile(out,'panel_ready_blind_summary.csv'));
fid=fopen(fullfile(out,'PLOT_ADAPTER_README.md'),'w'); fprintf(fid,'# Blind figure adapter\n\nLossless mapping of formal blind outputs. NaN denotes an unavailable quantity; no phase/noise-grid values are imputed. Legacy plotting scripts require a separate schema adapter after the full grid is regenerated.\n'); fclose(fid);
end
function t=mk(panel,protocol,id,et,eh,gt,gh,at,ah,rmse,status)
v={panel,protocol,string(id),string(et),string(eh),"NaN","NaN",string(gt),string(gh),string(at),string(ah),"NaN","NaN",string(rmse),string(status)}; t=cell2table(v,'VariableNames',{'figure_panel','protocol','case_id','eo_true','eo_hat','frequency_true_hz','frequency_hat_hz','gap_true_mm','gap_hat_mm','amplitude_true_mm','amplitude_hat_mm','phase_true_rad','phase_hat_rad','rmse_V','status'});
end
