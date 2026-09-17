function T = Run_03A_GlobalGapBudgetAudit()
%RUN_03A_GLOBALGAPBUDGETAUDIT Small audit before changing formal R3 runs.

thisDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(thisDir);
addpath(thisDir,'-begin');
addpath(fullfile(thisDir,'local_func'),'-begin');
R3 = R3_Protocol();
R3 = r3_build_context(mainDir,R3);
C0 = r3_freeze_reference_calibration(R3);

gaps = R3.P.gOperatingMm;
amps = [0 0.10 0.25 0.37];
repIds = 1:3;
rows = repmat(struct(),0,1);
for ig = 1:numel(gaps)
    for ia = 1:numel(amps)
        for ir = repIds
            seed = 2026087000 + 100*ig + 10*ia + ir;
            highMap = r3_generate_high_case(R3,C0,gaps(ig),amps(ia),seed);
            fits = r3_run_three_methods_global_audit(R3,C0,highMap,gaps(ig));
            row = r3_pack_result(R3,C0,highMap,fits,gaps(ig),amps(ia),ir,seed);
            row.audit_search_bounds_mm = [0.4 1.3];
            rows(end+1,1) = row; %#ok<AGROW>
            fprintf('global audit: g=%.2f A=%.2f rep=%d adaptive=%.2fs gHat=%.4f EO=%g\n', ...
                gaps(ig),amps(ia),ir,fits.adaptive.elapsed_r3_s,fits.adaptive.gHat,fits.adaptive.EO);
        end
    end
end
T = struct2table(rows);
outDir = fullfile(thisDir,'output','03a_global_gap_budget_audit');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T,fullfile(outDir,'global_gap_budget_audit.csv'));
save(fullfile(outDir,'global_gap_budget_audit.mat'),'T','R3','C0','-v7.3');
end
