function manifest = Run_EvidenceChain_20250527_20251222(varargin)
%RUN_EVIDENCECHAIN_20250527_20251222 Unified evidence-figure dispatcher.
%
% The dispatcher is unified, but the two experiments remain isolated: every
% case keeps its own result files, static bundle and output directory.  No
% identification is performed here; this entry point only reads frozen
% results and calls the existing diagnostic plotters.

p = inputParser;
addParameter(p,'Root',fileparts(mfilename('fullpath')),@ischar);
addParameter(p,'Cases',{'20250527','20251222'},@iscell);
addParameter(p,'Windows',1,@(x)isnumeric(x)&&isvector(x));
addParameter(p,'IncludeB5',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'MakeFigures',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});

root = p.Results.Root;
diagDir = fullfile(root,'diagnostics');
addpath(diagDir);
manifest = struct('schema','EVIDENCE_CHAIN_DISPATCH_V1', ...
    'created',char(datetime('now','Format','yyyyMMdd''T''HHmmss')),'root',root,'cases',[]);

for ic = 1:numel(p.Results.Cases)
    name = char(p.Results.Cases{ic});
    C = localCase(root,name,p.Results.IncludeB5);
    for ir = 1:numel(C.resultFiles)
      resultFile = C.resultFiles{ir}; runName = C.labels{ir};
      out = fullfile(root,'results',['evidence_chain_' name '_' runName]);
      if ~exist(out,'dir'), mkdir(out); end
      rec = struct('case',[name '_' runName],'resultFile',resultFile,'outputDir',out, ...
          'status','missing_result','figures',{{}},'errors',{{}});
      if ~isfile(resultFile)
        rec.errors = {sprintf('Result file not found: %s',resultFile)};
        manifest.cases = [manifest.cases; rec]; %#ok<AGROW>
        continue
      end
    rec.status = 'result_found';
    if p.Results.MakeFigures
        for iw = p.Results.Windows(:).'
            wout = fullfile(out,sprintf('window%02d',iw));
            if ~exist(wout,'dir'), mkdir(wout); end
            [rec,~] = localPlot(rec,@()Plot_B_DynamicSingleWindowAudit(resultFile,iw,wout,C.bundleFile), ...
                sprintf('B0_B2_window%02d',iw));
            [rec,~] = localPlot(rec,@()Plot_B1_StaticDynamicDecomposition(resultFile,iw,wout,C.bundleFile), ...
                sprintf('B1_decomposition_window%02d',iw));
            [rec,~] = localPlot(rec,@()Plot_B4_MethodAblationWaveform(resultFile,iw,wout,C.bundleFile), ...
                sprintf('B4_ablation_window%02d',iw));
        end
        % Batch summary is case-specific and therefore never combines the
        % 20250527 and 20251222 result bundles.
        if localHasResultStruct(resultFile)
            [rec,~] = localPlot(rec,@()Plot_C_DynamicBatchSummary(resultFile,fullfile(out,'batch')), 'C1_batch');
        else
            rec.errors{end+1} = 'C1_batch skipped: row-only result requires a case-specific comparison table.';
        end
    end
    rec.status = 'completed';
      manifest.cases = [manifest.cases; rec]; %#ok<AGROW>
    end
end

save(fullfile(root,'results','EvidenceChain_20250527_20251222_manifest.mat'),'manifest','-v7.3');
localWriteManifest(manifest,fullfile(root,'results','EvidenceChain_20250527_20251222_manifest.txt'));
end

function C = localCase(root,name,includeB5)
switch name
    case '20250527'
        base = fullfile(root,'20250527_gap_vibration_identification_v1','latest_programs_20260907','results');
        C.resultFiles = {fullfile(base,'r5_sensor_conditioned_dynamic_fullwave_20260914.mat')};
        C.labels = {'R5'};
        C.bundleFile = fullfile(base,'fixed_gap','20250526_2500-3500_t400','FixedGap_B1_S136_T001p5s.mat');
    case '20251222'
        base = fullfile(root,'20251222_gap_vibration_identification_v1','latest_programs_20260905','results','gap_aware');
        % B1 and B5 are intentionally separate experiments; B1 is the
        % default entry and B5 can be dispatched independently by rerunning
        % with a case-specific result file in this adapter.
        C.resultFiles = {fullfile(base,'Main_GapAware_VPTopK_FullWave_20251222_B1_S123_r01_gapaware.mat')};
        C.labels = {'B1_R01'};
        C.bundleFile = '';
        if includeB5
            C.resultFiles{end+1} = fullfile(base,'Main_GapAware_VPTopK_FullWave_20251222_B5_S123_r04_gapaware.mat');
            C.labels{end+1} = 'B5_R04';
        end
    otherwise
        error('EvidenceChain:UnknownCase','Unknown case %s.',name);
end
end

function tf = localHasResultStruct(file)
% Detect the public Result schema without loading large waveform arrays.
wh = whos('-file',file); names = {wh.name}; tf = any(strcmp(names,'Result'));
end

function [rec,ok] = localPlot(rec,fun,label)
ok = false;
try
    fun(); rec.figures{end+1} = label; ok = true;
catch ME
    rec.errors{end+1} = sprintf('%s: %s',label,ME.message);
end
end

function localWriteManifest(M,file)
fid=fopen(file,'w'); if fid<0, return; end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'schema=%s\ncreated=%s\n',M.schema,M.created);
for k=1:numel(M.cases)
    c=M.cases(k); fprintf(fid,'case=%s\tstatus=%s\tresult=%s\n',c.case,c.status,c.resultFile);
    for j=1:numel(c.figures), fprintf(fid,'  figure=%s\n',c.figures{j}); end
    for j=1:numel(c.errors), fprintf(fid,'  error=%s\n',c.errors{j}); end
end
end
