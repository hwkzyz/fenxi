%% Run_NewFlow_20241106
% One-command runner for the current numbering-first workflow.
% Step00-Step05 tune parameters at the top of each script.
%
% Current main route:
%   Step05 = low-speed template calibration
%   Step06 = direct-template identification with inline high-speed extraction
%   Step07 = audit
%
% Visualization / scan / deep-audit helpers were moved into:
%   diagnostics/

clear; close all; clc;

routeDir = fileparts(mfilename('fullpath'));
cd(routeDir);

Step00_BuildLowSpeedSensorConfig_20241106;
Step01_BuildLowSpeedReferenceFingerprint_20241106;
Step02_NumberLowSpeedSensors_20241106;
Step03_SelectHighSpeedRegion_20241106;
Step04P_BuildHighSpeedPeakCache_20241106;
Step04_NumberHighSpeedSensorsInRegion_20241106;
Step05_BuildLowSpeedTemplateLibrary_20241106;
Step06_RunIdentificationByBlade_20241106;
Step07_AuditNumberingAndCorrespondence_20241106;
