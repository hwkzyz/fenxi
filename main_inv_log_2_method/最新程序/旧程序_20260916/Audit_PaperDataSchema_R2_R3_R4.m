function Report = Audit_PaperDataSchema_R2_R3_R4()
%AUDIT_PAPERDATASCHEMA_R2_R3_R4 Read-only audit before paper-data replacement.
% No figure, source-data, or formal result file is modified.
root=fileparts(mfilename('fullpath'));
files={ ...
 fullfile(root,'r2_profiled_recovery','output','unified_blind_backend','r2_unified_blind_backend.csv'), ...
 fullfile(root,'r3_heldout_recovery','output','unified_blind_backend','r3_unified_blind_backend.csv'), ...
 fullfile(root,'r4_pc1_manifold_transfer','blind_unknown_gap','output','gaponly_vp_unknown_eo','case_results.csv'), ...
 fullfile(root,'r3_heldout_recovery','output','05_main_15db','main_detail_merged.csv'), ...
 fullfile(root,'r4_pc1_manifold_transfer','output','transferred_surface_dynamic_validation','transferred_surface_dynamic_results.csv')};
names={"R2_unified";"R3_unified";"R4_unified";"R3_legacy_detail";"R4_legacy_dynamic"};
rows=repmat(struct('dataset',"",'exists',false,'rows',NaN,'columns',NaN,'hasTruthEO',false,'hasEstimatedEO',false,'hasGapEstimate',false,'hasAmplitudeEstimate',false,'hasIdentifiability',false,'note',""),numel(files),1);
for i=1:numel(files)
 rows(i).dataset=names{i}; rows(i).exists=isfile(files{i});
 if ~rows(i).exists, rows(i).note="missing"; continue; end
 T=readtable(files{i},'TextType','string'); v=string(T.Properties.VariableNames); rows(i).rows=height(T); rows(i).columns=width(T);
 rows(i).hasTruthEO=any(v=="eo_true"); rows(i).hasEstimatedEO=any(v=="eo_hat"|v=="eo_est");
 rows(i).hasGapEstimate=any(v=="gap_hat_mm"|v=="g_est_mm"); rows(i).hasAmplitudeEstimate=any(v=="A_hat_mm"|v=="A_est_mm"); rows(i).hasIdentifiability=any(v=="identifiability_status");
 rows(i).note=string(files{i});
end
Report=struct2table(rows); out=fullfile(root,'paper_figures','updated_source_data'); if ~exist(out,'dir'),mkdir(out);end
writetable(Report,fullfile(out,'paper_data_schema_audit.csv')); save(fullfile(out,'paper_data_schema_audit.mat'),'Report'); disp(Report)
end
