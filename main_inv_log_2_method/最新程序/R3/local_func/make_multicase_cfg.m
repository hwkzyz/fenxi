function cfgAna = make_multicase_cfg(cfg)
%make_multicase_cfg  Match the saved Step 6G multi-scenario settings.

cfgAna = cfg;
cfgAna.fs = 2e6;
cfgAna.NumRevs_high = 4;
cfgAna.gapQueryN = 241;
cfgAna.dxStaticGrid = linspace(-0.6, 0.6, 161);
cfgAna.f1Grid = 400:5:650;
cfgAna.f2Grid = 1100:5:1400;
% Physical frequency metadata for structured S/A/SS/SA analyses. The
% legacy split grids above remain unchanged for production compatibility.
cfgAna = define_frequency_range(cfgAna, [100, 1600], 2);
cfgAna.singleFreqGrid = cfgAna.frequencyGridFromRangeHz(:).';
cfgAna.singleAsyncTopK = 3;
cfgAna.singleAsyncMaxTopK = 5;
cfgAna.singleAsyncLowSNRThreshold = 20;
cfgAna.singleAsyncLowSNRTopK = 5;
cfgAna.singleAsyncVeryLowSNRThreshold = 10;
cfgAna.singleAsyncVeryLowSNRTopK = 8;
cfgAna.singleAsyncExpandRhoThreshold = 0.05;
cfgAna.singleAsyncCandidateMinSeparationHz = 20;
cfgAna.singleAsyncRefineHalfWidthHz = 10;
cfgAna.singleAsyncFineStepHz = 0.1;
% Keep two asynchronous local candidates before the complete voltage-domain
% fit.  A single inverse-map winner is not reliable at low SNR; the final
% forward residual still selects the production result.
cfgAna.singleSAFineKeep = 2;
cfgAna.singleSALowSNRThreshold = 20;
cfgAna.singleSALowSNRTopK = 5;
cfgAna.singleSALowSNRFineKeep = 3;
cfgAna.singleSAVeryLowSNRThreshold = 10;
cfgAna.singleSAVeryLowSNRTopK = 8;
cfgAna.singleSAVeryLowSNRFineKeep = 4;
% The inverse interval must cover the formal 0.5 mm vibration range plus
% the allowed 0.2 mm offset. It is only a candidate generator; the final
% voltage fit still uses the wider 1.5 mm numerical amplitude bound.
cfgAna.inverseMapHalfWidthMm = 0.70;
cfgAna.numVarproCandidates = 10;
cfgAna.gapStaticTargetSamples = 4000;
cfgAna.gapStaticCoarseGapN = 41;
cfgAna.gapStaticCoarseDxN = 41;
cfgAna.gapStaticFineGapN = 31;
cfgAna.gapStaticFineDxN = 31;
cfgAna.gapStaticFineGapHalfWidth = 0.03;
cfgAna.gapStaticFineDxHalfWidth = 0.06;
cfgAna.gapProjTargetSamples = 5000;
cfgAna.gapProjCoarseN = 61;
cfgAna.gapProjFineN = 31;
cfgAna.gapProjHalfWidth = 0.08;
cfgAna.gapProjFineHalfWidth = 0.02;
cfgAna.gapProjectionLambda = 1;
cfgAna.gapInitAutoFreqCount = 4;
cfgAna.gapInitAutoFreqMinSep = 80;
cfgAna.confGapLocalHalfWidth = 0.03;
cfgAna.confGapConsensusN = 11;
cfgAna.confGapExpandTol = 0.02;
cfgAna.confGapCoarseN = 7;
cfgAna.confTopGapKeep = 3;
cfgAna.confFreqVoteTolHz = 8;
cfgAna.confRhoJWeight = 0.5;
cfgAna.confUseProjectedJointMuList = [0, 0.15, 0.75];
cfgAna.confRefineGapHalfWindow = 1;
cfgAna.confRefineGapKeep = 3;
cfgAna.mainVpGapHalfWidth = 0.06;
cfgAna.mainVpGapN = 61;
cfgAna.mainVpCoarseCount = 21;
cfgAna.mainVpKeepCoarse = 5;
cfgAna.mainVpRelTol = 0.02;
cfgAna.mainVpMaxKeep = 7;
cfgAna.mainVpFineHalfWindow = 5;
cfgAna.mainVpRefineCount = 5;
cfgAna.mainUseProjectedJointMuList = [0, 0.15, 0.75];
end
