function M = Refresh_Fig02_Fig03_BlindSourceData()
%REFRESH_FIG02_FIG03_BLINDSOURCEDATA Stage new blind data for figure review.
% This routine never edits the legacy fig02/fig03 source_data folders.
root = fileparts(mfilename('fullpath'));
batch = fullfile(root,'paper_figures','updated_source_data','formal_blind_batch');
assert(isfolder(batch),'Run Run_Formal_Blind_R2_R3_R4_Batch first.');
dest = fullfile(root,'paper_figures','updated_source_data','fig02_fig03_blind_review');
if ~exist(dest,'dir'), mkdir(dest); end
files = {fullfile(batch,'R2_formal_blind.csv'), ...
         fullfile(batch,'R3_formal_blind.csv'), ...
         fullfile(batch,'R4_formal_blind.csv')};
names = {'Fig02_R2_blind_joint_recovery.csv', ...
         'Fig02_R3_blind_near_synchronous.csv', ...
         'Fig03_R4_blind_dynamic_geometry.csv'};
for k=1:numel(files)
    if isfile(files{k}), copyfile(files{k},fullfile(dest,names{k}),'f'); end
end
manifest = table(string(names(:)), isfile(files(:)), ...
    'VariableNames',{'file','available'});
writetable(manifest,fullfile(dest,'manifest.csv'));
fid=fopen(fullfile(dest,'README.md'),'w');
fprintf(fid,'# Updated Fig. 2 / Fig. 3 blind-review data\n\n');
fprintf(fid,'Generated without modifying legacy figures or source data.\n');
fprintf(fid,'EO estimation scans candidates 1:30; truth parameters are not estimator inputs.\n');
fprintf(fid,'Fig. 2 panels d--i use R2/R3 files; Fig. 3f uses the R4 file.\n');
fprintf(fid,'The existing plotting scripts still require a schema adapter before plotting these files.\n');
fclose(fid);
M=struct('directory',dest,'manifest',manifest);
end
