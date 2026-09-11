%% Independent static waveform-manifold analysis runner
clear; clc; close all;
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
run('Step01_BuildStaticWaveformFamily.m');
run('Step02_BuildBladeWaveformSurfaces.m');
run('Step06_V1_SameBladeTransition_20251222.m');
run('Step07_V2_LeaveOneGap_20251222.m');
run('Step03_EstimateLocalGapPaths_20251222.m');
run('Step04_LeaveOneBladeLocalization_20251222.m');
run('Step05_PeakShiftDiagnostics_20251222.m');
run('Step09_V4_ExperimentalClosure_20251222.m');
run('Step10_ValidationGate_20251222.m');
run('Step11_CompareMethods_V1V4_20251222.m');
run('Step12_CompareIncrementRoutes_V1V3_20251222.m');
run('Step13_V3_ManifoldLocalization_20251222.m');
run('Step13_VisualizeValidationFourLevels_20251222.m');
run('Step14_PlotValidationWaveformComparisons_20251222.m');
