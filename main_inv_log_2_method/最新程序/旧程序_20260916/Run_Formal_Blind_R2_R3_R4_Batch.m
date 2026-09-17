function S = Run_Formal_Blind_R2_R3_R4_Batch(opts)
%RUN_FORMAL_BLIND_R2_R3_R4_BATCH Re-run the frozen blind protocols.
% Truth parameters are consumed only by each waveform generator; all three
% estimators perform their own EO search.  Results are written to a new
% paper_figures/updated_source_data/formal_blind_batch folder.
if nargin < 1, opts = struct(); end
root = fileparts(mfilename('fullpath'));
if ~isfield(opts,'r2Cases'), opts.r2Cases = 6; end
if ~isfield(opts,'r3Cases'), opts.r3Cases = 3; end
if ~isfield(opts,'runR4'), opts.runR4 = true; end
out = fullfile(root,'paper_figures','updated_source_data','formal_blind_batch');
if ~exist(out,'dir'), mkdir(out); end
assert(exist(fullfile(root,'local_func','run_unified_blind_dynamic_backend.m'),'file')==2, ...
    'Unified blind backend is missing from the latest program.');

% Register only formal folders; the entry-point guards reject shadowed code.
addpath(root,'-begin');
R2 = Run_R2_BlindSingleFrequency(opts.r2Cases);
R3 = Run_R3_BlindSingleFrequency(opts.r3Cases);
assert(all(R2.table.eo_true>=1 & R2.table.eo_true<=30 & R2.table.eo_hat>=1 & R2.table.eo_hat<=30));
assert(all(R3.table.eo_true>=1 & R3.table.eo_true<=30 & R3.table.eo_hat>=1 & R3.table.eo_hat<=30));
writetable(R2.table,fullfile(out,'R2_formal_blind.csv'));
writetable(R3.table,fullfile(out,'R3_formal_blind.csv'));

R4 = struct();
if opts.runR4
    r4out = fullfile(out,'R4');
    o = struct('outputDir',r4out);
    R4 = Main_10_UnifiedBlindBackendValidation(o);
    assert(all(strlength(string(R4.cases.identification_method))>0));
    writetable(R4.cases,fullfile(out,'R4_formal_blind.csv'));
end

summary = struct();
summary.R2_eo_recovery = mean(R2.table.eo_true == R2.table.eo_hat);
summary.R2_identified = mean(string(R2.table.status)=="identified");
summary.R3_eo_recovery = mean(R3.table.eo_true == R3.table.eo_hat);
summary.R3_identified = mean(string(R3.table.status)=="identified");
if opts.runR4
    summary.R4_eo_recovery = mean(string(R4.cases.eo_true)==string(R4.cases.eo_hat));
    summary.R4_identified = mean(string(R4.cases.identifiability_status)=="identified");
else
    summary.R4_eo_recovery = NaN; summary.R4_identified = NaN;
end
summary.generated_at = string(datetime('now'));
S = struct('R2',R2,'R3',R3,'R4',R4,'summary',summary,'outputDir',out);
save(fullfile(out,'formal_blind_batch.mat'),'S','-v7.3');
fid = fopen(fullfile(out,'formal_blind_summary.md'),'w');
fprintf(fid,'# Formal blind R2--R4 batch\n\n');
fprintf(fid,'All EO candidates were searched blindly over 1:30.\n\n');
fprintf(fid,'Truth parameters are isolated to waveform generation and final error evaluation; no truth EO/gap field is passed to the estimator.\n\n');
fprintf(fid,'| Protocol | EO recovery | identified |\n|---|---:|---:|\n');
fprintf(fid,'| R2 | %.3f | %.3f |\n',summary.R2_eo_recovery,summary.R2_identified);
fprintf(fid,'| R3 | %.3f | %.3f |\n',summary.R3_eo_recovery,summary.R3_identified);
if opts.runR4, fprintf(fid,'| R4 | %.3f | %.3f |\n',summary.R4_eo_recovery,summary.R4_identified); end
fclose(fid);
disp(summary);
end
