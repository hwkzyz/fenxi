function Result = Merge_05_Main15dB_Workers()
%MERGE_05_MAIN15DB_WORKERS Audit and merge the formal and worker checkpoints.

root=fileparts(mfilename('fullpath'));P=R3_Protocol();outDir=fullfile(root,'output','05_main_15db');
files=[string(fullfile(outDir,'main_detail_checkpoint.csv'));...
    string(fullfile(root,'output','05_main_15db_workers','worker_01','detail_checkpoint.csv'));...
    string(fullfile(root,'output','05_main_15db_workers','worker_02','detail_checkpoint.csv'));...
    string(fullfile(root,'output','05_main_15db_workers','worker_03','detail_checkpoint.csv'))];
if any(~isfile(files)),error('r3:MissingWorkerOutput','All formal and worker checkpoints are required before merging.');end
parts=cell(numel(files),1);for i=1:numel(files),parts{i}=readtable(files(i),'TextType','string');end
Detail=vertcat(parts{:});audit_complete_matrix(Detail,P);
[Summary,Paired]=Summarize_R3(Detail,P);
writetable(Detail,fullfile(outDir,'main_detail_merged.csv'));
writetable(Summary,fullfile(outDir,'main_summary.csv'));
writetable(Paired,fullfile(outDir,'main_paired_summary.csv'));
Result=struct('P',P,'Detail',Detail,'Summary',Summary,'Paired',Paired,'source_files',files,'outputDir',outDir);
save(fullfile(outDir,'main_result.mat'),'Result','-v7.3');
end

function audit_complete_matrix(T,P)
required=3*P.replicateCountMaximum*numel(P.gMainMm)*(1+numel(P.amplitudeMainMm));
if height(T)~=required,error('r3:RowCount','Expected %d rows, found %d.',required,height(T));end
key=T(:,{'method','g_truth_mm','A_true_mm','high_seed'});
if height(unique(key,'rows'))~=height(T),error('r3:DuplicateRows','Duplicate method/gap/amplitude/seed rows found.');end
for g=P.gMainMm
    for A=[P.amplitudeNullMm P.amplitudeMainMm]
        q=abs(T.g_truth_mm-g)<1e-12&T.A_true_mm==A;
        if nnz(q)~=3*P.replicateCountMaximum,error('r3:IncompleteCell','Cell g=%.3f, A=%.3f has %d rows.',g,A,nnz(q));end
        for seed=P.mainHighSeeds
            if nnz(q&T.high_seed==seed)~=3,error('r3:IncompleteTriplet','Cell g=%.3f, A=%.3f has an incomplete method triplet.',g,A);end
        end
    end
end
names=["fixed","adaptive","state_matched"];
if ~isequal(sort(unique(T.method)),sort(names.')),error('r3:MethodSet','Unexpected method set.');end
end
