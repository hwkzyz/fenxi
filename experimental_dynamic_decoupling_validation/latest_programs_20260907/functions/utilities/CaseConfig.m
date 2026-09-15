function C = CaseConfig()
%CASECONFIG Compatibility view of the frozen package configuration.
P = ProjectionFlow_Config_20250527();
C = struct();
C.dataset = P.dataset;
C.bladeId = P.identification.targetBlade;
C.sensorIds = P.identification.analysisSensors;
C.gapSensors = P.identification.gapSensors;
C.sensorTag = ['S', sprintf('%d', C.sensorIds)];
C.caseTag = sprintf('B%d_%s', C.bladeId, C.sensorTag);
C.flowConfig = P;
end
