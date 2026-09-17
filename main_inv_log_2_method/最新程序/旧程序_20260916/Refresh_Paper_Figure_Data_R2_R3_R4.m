function Manifest = Refresh_Paper_Figure_Data_R2_R3_R4(doRerun)
%REFRESH_PAPER_FIGURE_DATA_R2_R3_R4
% Re-run the unified blind R2/R3/R4 simulations and publish their results
% beside the existing paper data without changing any plotting program.
%
% Existing formal paper MAT files are not overwritten when their schema is
% tied to the historical fixed-EO experiment.  A manifest records that
% distinction so that a figure cannot silently mix blind and oracle data.

here = fileparts(mfilename('fullpath'));
if nargin<1, doRerun=true; end
oldRoot = fullfile(here,'..','旧版程序');
stamp = datestr(now,'yyyymmdd_HHMMSS');
backupRoot = fullfile(oldRoot,['paper_data_backup_' stamp]);
if ~exist(backupRoot,'dir'), mkdir(backupRoot); end

% Re-run the three canonical blind simulations.
addpath(here,'-begin');
addpath(fullfile(here,'r4_pc1_manifold_transfer','blind_unknown_gap'),'-begin');
if doRerun
    R2 = Run_R2_BlindSingleFrequency(6);
    R3 = Run_R3_BlindSingleFrequency(3);
    R4 = Main_09_BlindGapOnlyVPDynamicValidation();
else
    R2 = struct('outputDir',fullfile(here,'r2_profiled_recovery','output','unified_blind_backend'));
    R3 = struct('outputDir',fullfile(here,'r3_heldout_recovery','output','unified_blind_backend'));
    R4 = struct('output_dir',fullfile(here,'r4_pc1_manifold_transfer','blind_unknown_gap','output','gaponly_vp_unknown_eo'));
end

% Keep a paper-facing, non-destructive data package.
paperRoot = fullfile(here,'paper_figures');
dataRoot = fullfile(paperRoot,'updated_source_data');
if ~exist(dataRoot,'dir'), mkdir(dataRoot); end
copyfile(fullfile(R2.outputDir,'r2_unified_blind_backend.csv'), ...
    fullfile(dataRoot,'R2_unified_blind_backend.csv'),'f');
copyfile(fullfile(R3.outputDir,'r3_unified_blind_backend.csv'), ...
    fullfile(dataRoot,'R3_unified_blind_backend.csv'),'f');
copyfile(fullfile(char(R4.output_dir),'case_results.csv'), ...
    fullfile(dataRoot,'R4_unified_blind_backend.csv'),'f');

% Freeze the R2 source files consumed by Plot_Fig02 (schema-compatible).
fig2 = fullfile(paperRoot,'fig02');
if exist(fullfile(fig2,'Freeze_Fig02_SourceData.m'),'file')
    old = pwd; c = onCleanup(@()cd(old)); cd(fig2);
    Freeze_Fig02_SourceData(); clear c
end

% Back up historical formal files before any future manual synchronization.
formalFiles = { ...
    fullfile(here,'r3_heldout_recovery','output','18_paper_data_sync','r3_paper_result.mat'), ...
    fullfile(here,'r3_heldout_recovery','output','18_paper_data_sync','r3_paper_summary.csv'), ...
    fullfile(here,'r4_pc1_manifold_transfer','output','transferred_surface_dynamic_validation','transferred_surface_dynamic_results.csv')};
for i=1:numel(formalFiles)
    f = formalFiles{i};
    if isfile(f)
        [~,name,ext] = fileparts(f);
        copyfile(f,fullfile(backupRoot,[name ext]),'f');
    end
end

% Machine-readable contract for the unchanged plotting programs.
Manifest = table( ...
    string({'R2';'R3';'R4';'Fig02_static';'Fig03_formal';'Fig04_geometry'}), ...
    string({'blind_unified';'blind_unified';'blind_unified';'schema_compatible';'historical_schema';'historical_schema'}), ...
    string({fullfile(dataRoot,'R2_unified_blind_backend.csv'); ...
            fullfile(dataRoot,'R3_unified_blind_backend.csv'); ...
            fullfile(dataRoot,'R4_unified_blind_backend.csv'); ...
            fullfile(fig2,'source_data','fig02_bc_canonical_profile.csv'); ...
            fullfile(here,'r3_heldout_recovery','output','18_paper_data_sync','r3_paper_result.mat'); ...
            fullfile(here,'r4_pc1_manifold_transfer','output','transferred_surface_dynamic_validation','transferred_surface_dynamic_results.csv')}), ...
    string({'latest blind rerun';'latest blind rerun';'latest blind rerun';'static profile retained';'not replaced: fixed-EO schema';'not replaced: geometry-transfer schema'}), ...
    'VariableNames',{'dataset','contract','path','note'});
writetable(Manifest,fullfile(dataRoot,'paper_data_refresh_manifest.csv'));
save(fullfile(dataRoot,'paper_data_refresh_manifest.mat'),'Manifest','-v7.3');
fprintf('PAPER_DATA_REFRESH_OK backup=%s data=%s\n',backupRoot,dataRoot);
end
