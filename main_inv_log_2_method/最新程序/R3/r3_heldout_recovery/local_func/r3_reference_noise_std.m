function [sigmaV,info] = r3_reference_noise_std(R3,snrDb)
%R3_REFERENCE_NOISE_STD Freeze voltage noise from the 0.10-mm reference motion.

P=R3.P;cfg=R3.cfg;Fraw=r3_raw_response(R3,P.gReferenceMm);
fTrue=P.eoTrue*cfg.RPM_high/60;
vib=simulate_rotating_waveform_from_template(Fraw,cfg.RPM_high,...
    cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,R3.libCal.domain,...
    P.referenceAmplitudeMm,fTrue,P.referencePhaseRad,'fixed_std');
stat=simulate_rotating_waveform_from_template(Fraw,cfg.RPM_high,...
    cfg.NumRevs_high,cfg.fs,cfg.R_tip,cfg.alpha_k,0,R3.libCal.domain,...
    [],[],[],'fixed_std');
vib.V_cap=vib.V_clean;stat.V_cap=stat.V_clean;
mv=map_highspeed_to_space(vib,cfg.alpha_k,cfg.R_tip,R3.libCal.domain,.02);
ms=map_highspeed_to_space(stat,cfg.alpha_k,cfg.R_tip,R3.libCal.domain,.02);
if numel(mv.V_a)~=numel(ms.V_a)
    error('r3:NoiseReferenceAlignment','Reference dynamic/static maps are misaligned.');
end
referenceVoltageRms=rms(mv.V_a-ms.V_a);
sigmaV=referenceVoltageRms/10^(snrDb/20);
info=struct('snr_db',snrDb,'reference_amplitude_mm',P.referenceAmplitudeMm,...
    'reference_gap_mm',P.gReferenceMm,'reference_voltage_rms_V',referenceVoltageRms,...
    'noise_std_V',sigmaV,'definition',P.referenceSnrDefinition);
end
