function Result=Run_04_X01_Evidence_Synthesis()
%RUN_04_X01_EVIDENCE_SYNTHESIS Import frozen X01 evidence without rerunning it.
analysisDir=fileparts(mfilename('fullpath'));programDir=fileparts(analysisDir);
srcDir=fullfile(programDir,'low_high_gap_fast_bridge','output');
files={fullfile(srcDir,'x01_deterministic_bias_audit','summary.csv'), ...
       fullfile(srcDir,'x01_template_stage_ablation','summary.csv'), ...
       fullfile(srcDir,'x01_spatial_discretization_audit','summary.csv')};
exists=cellfun(@isfile,files);
if ~exists(1),error('identifiability:MissingX01','Frozen X01 summary is missing.');end
Audit=readtable(files{1},'TextType','string');
formal=Audit(Audit.model=="formal_adaptive_template",:);
exact=Audit(Audit.model=="exact_low_direct_static",:);
Summary=table();
Summary.formal_wrong_minus_true_rmse_V=formal.wrong_minus_true_oracle_rmse_V;
Summary.formal_min_angle_deg=formal.min_angle_deg;
Summary.formal_rho_max=formal.rho_max;
Summary.formal_true_parameter_rmse_V=formal.true_parameter_rmse_V;
Summary.exact_wrong_minus_true_rmse_V=exact.wrong_minus_true_oracle_rmse_V;
Summary.exact_true_parameter_rmse_V=exact.true_parameter_rmse_V;
Summary.bias_flip_present=Summary.formal_wrong_minus_true_rmse_V<0 & ...
    Summary.exact_wrong_minus_true_rmse_V>0;
Summary.near_collinear=Summary.formal_rho_max>.99;
Summary.supports_model_separation_bias_competition= ...
    Summary.bias_flip_present & Summary.near_collinear;
out=fullfile(analysisDir,'output','04_x01_evidence_synthesis');
if ~exist(out,'dir'),mkdir(out);end
writetable(Audit,fullfile(out,'imported_x01_audit.csv'));
writetable(Summary,fullfile(out,'evidence_summary.csv'));
Manifest=table(string(files(:)),exists(:),'VariableNames',{'source_file','exists'});
writetable(Manifest,fullfile(out,'source_manifest.csv'));
write_report(out,Summary,Manifest);
disp(Summary);Result=struct('Audit',Audit,'Summary',Summary,'Manifest',Manifest,'outputDir',out);
end
function write_report(out,S,M)
fid=fopen(fullfile(out,'README.md'),'w');fprintf(fid,'# X01 evidence synthesis\n\n');
fprintf(fid,'This report imports frozen evidence and does not rerun or tune the solver.\n\n');
fprintf(fid,'- Formal wrong-minus-true oracle RMSE: %.9g V\n',S.formal_wrong_minus_true_rmse_V);
fprintf(fid,'- Exact-model wrong-minus-true oracle RMSE: %.9g V\n',S.exact_wrong_minus_true_rmse_V);
fprintf(fid,'- Conditional canonical correlation: %.9g\n',S.formal_rho_max);
fprintf(fid,'- Minimum principal angle: %.9g deg\n',S.formal_min_angle_deg);
fprintf(fid,'- Supports separation-bias competition for X01: %d\n\n',S.supports_model_separation_bias_competition);
fprintf(fid,'This is direct evidence for X01 only. General prediction requires locked metrics on new EO pairs.\n\n');
fprintf(fid,'Available source artifacts: %d/%d.\n',sum(M.exists),height(M));fclose(fid);
end
