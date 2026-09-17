function T = Run_03A_GlobalGapAudit(replicateIds)
%RUN_03A_GLOBALGAPAUDIT Validate the global adaptive search before formal R3.

if nargin < 1, replicateIds = 1; end
thisDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(thisDir);
addpath(thisDir,'-begin');
addpath(fullfile(thisDir,'local_func'),'-begin');
R3 = r3_build_context(mainDir,R3_Protocol());
[sigmaV,~] = r3_reference_noise_std(R3,R3.P.referenceSnrDb);
[C0,~] = r3_freeze_reference_calibration(R3,sigmaV,R3.P.mainLowSeed);

gaps = R3.P.gOperatingMm;
amps = [0.10 0.37];
records = repmat(struct('method',"",'gTruthMm',NaN,'AtrueMm',NaN,'rep',NaN,...
    'gEstMm',NaN,'AEstMm',NaN,'EO',NaN,'phaseErrorRad',NaN,...
    'elapsedS',NaN,'vibSuccess',false),0,1);
for ig = 1:numel(gaps)
    for ia = 1:numel(amps)
        for ir = replicateIds
            seed = R3.P.mainHighSeeds(ir);
            phi = R3.P.phaseValuesRad(ir);
            [~,highMap,truth] = r3_generate_high_case(R3,gaps(ig),amps(ia),phi,sigmaV,seed);
            fits = r3_run_three_methods_global_final(R3,C0,highMap,gaps(ig));
            names = {'fixed','adaptive','state_matched'};
            for im = 1:numel(names)
                fit = fits.(names{im});
                phaseError = abs(atan2(sin(fit.phi_id-truth.phi_mapped_rad),...
                    cos(fit.phi_id-truth.phi_mapped_rad)));
                ampTol = max(R3.P.amplitudeAbsToleranceMm,...
                    R3.P.amplitudeRelTolerance*truth.A_truth_mapped_mm);
                success = fit.eo_id == truth.eo_true && ...
                    abs(fit.A_id-truth.A_truth_mapped_mm) <= ampTol && ...
                    phaseError <= R3.P.phaseToleranceRad;
                records(end+1) = struct('method',string(names{im}),... %#ok<AGROW>
                    'gTruthMm',gaps(ig),'AtrueMm',amps(ia),'rep',ir,...
                    'gEstMm',fit.g_used,'AEstMm',fit.A_id,'EO',fit.eo_id,...
                    'phaseErrorRad',phaseError,'elapsedS',fit.elapsed_r3_s,...
                    'vibSuccess',success);
            end
            fprintf('g=%.2f A=%.2f rep=%d: adaptive g=%.4f EO=%g A=%.4f, %.2fs\n',...
                gaps(ig),amps(ia),ir,fits.adaptive.g_used,fits.adaptive.eo_id,...
                fits.adaptive.A_id,fits.adaptive.elapsed_r3_s);
        end
    end
end
T = struct2table(records);
outDir = fullfile(thisDir,'output','03a_global_gap_audit');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T,fullfile(outDir,'global_gap_audit.csv'));
save(fullfile(outDir,'global_gap_audit.mat'),'T','R3','C0','-v7.3');
end
