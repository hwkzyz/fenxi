%% Run_NewFlow_20251222
% One-command runner for the 20251222 NewFlow route.

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
cd(routeDir);

Step00_BuildLowSpeedSensorConfig_20251222;
Step01_BuildLowSpeedReferenceFingerprint_20251222;
Step02_NumberLowSpeedSensors_20251222;
Step03_SelectHighSpeedRegion_20251222;
Step04P_BuildHighSpeedPeakCache_20251222;
Step04_NumberHighSpeedSensorsInRegion_20251222;
Step05_BuildLowSpeedTemplateLibrary_20251222;
Step06_RunIdentificationByBlade_20251222;
Step07_AuditNumberingAndCorrespondence_20251222;
