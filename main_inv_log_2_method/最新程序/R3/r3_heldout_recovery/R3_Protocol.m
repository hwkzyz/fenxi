function P = R3_Protocol()
%R3_PROTOCOL Frozen protocol for the COMSOL-based held-out R3 analysis.

P.version = "R3-COMSOL-heldout-v3-global-gap-20260826";
P.claim = "Conditional on a frozen low-speed calibration, clearance-state adaptation prevents held-out operating-state mismatch from being attributed to synchronous blade motion.";

P.gAllMm = 0.2:0.1:1.5;
P.gCalMm = [0.2 0.3 0.4 0.6 0.8 1.1 1.3 1.4 1.5];
P.gOperatingMm = [0.5 0.7 0.9 1.0 1.2];
P.gReferenceMm = 0.8;
P.gMainMm = [P.gOperatingMm P.gReferenceMm];
P.gapBoundsMm = [0.2 1.5];
P.adaptiveSearchBoundsMm = [0.4 1.3];
P.adaptiveCoarseGapGridMm = 0.4:0.1:1.3;
P.adaptiveSearchScope = "tested numerical search range, not an engineering maximum";

P.amplitudeMainMm = [0.10 0.15 0.20 0.25 0.30 0.37];
P.amplitudeNullMm = 0;
P.eoTrue = 10;
P.deltaFTrueHz = 0;
P.frequencyRangeHz = [300 1500];
P.phaseValuesRad = 2*pi*(((1:75)-0.5)/75);

P.referenceSnrDb = 15;
P.supplementarySnrDb = [10 20];
P.referenceAmplitudeMm = 0.10;
P.referencePhaseRad = pi/4;
P.referenceSnrDefinition = "reference-amplitude-equivalent voltage SNR";

P.amplitudeAbsToleranceMm = 0.020;
P.amplitudeRelTolerance = 0.10;
P.phaseToleranceRad = deg2rad(20);
P.gapToleranceMm = 0.020;
P.targetProbability = 0.90;
P.wilsonAlpha = 0.05;

P.replicateCountInitial = 25;
P.replicateCountMaximum = 75;
P.initialReplicateIndices = 1:3:73;
P.remainingReplicateIndices = setdiff(1:75,P.initialReplicateIndices,'stable');
P.mainLowSeed = 2026082501;
P.mainHighSeeds = 2026082600 + (1:75);
P.detectionHighSeeds = 2026083000 + (1:500);
P.lowSensitivitySeeds = 2026084000 + (1:10);
P.detectionCalibrationCount = 500;
P.detectionQuantile = 0.95;

P.smokeGapsMm = [0.5 0.9 0.8 1.2];
P.smokeAmplitudesMm = [0 0.10 0.25 0.37];
P.smokeReplicateIndices = 1:3;

P.calibratedStateMatchedGateAmplitudeMm = [0 0.10 0.25 0.37];
P.calibratedStateMatchedMaxInterpNrmse = 0.05;

P.saveRepresentativeWaveforms = true;
end
