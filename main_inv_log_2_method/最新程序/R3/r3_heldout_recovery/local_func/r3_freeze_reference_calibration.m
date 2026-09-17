function [C0,low0] = r3_freeze_reference_calibration(R3,sigmaV,lowSeed)
%R3_FREEZE_REFERENCE_CALIBRATION Build the one predeclared main calibration.

Fref=r3_raw_response(R3,R3.P.gReferenceMm);
low0=simulate_low_speed_template(Fref,R3.cfg,R3.libCal.domain,...
    sigmaV,lowSeed,'fixed_std');
C0=calibrate_inv_log_2_low_speed(low0,R3.libCal,R3.cfg);
C0.r3_truth_reference_gap_mm=R3.P.gReferenceMm;
C0.r3_low_seed=lowSeed;
C0.r3_noise_std_V=sigmaV;
end
