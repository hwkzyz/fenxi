function Fits = r3_run_three_methods(R3,C0,highMap,gTruth)
%R3_RUN_THREE_METHODS Run paired fixed, adaptive and calibrated state-matched fits.

cfg=R3.cfg;
cfg.route30FrequencyStructureMode="single_sync";

ticFit=tic;
adaptive=run_single_near_sync_voltage_vp(highMap,C0.templateModel,cfg,C0.lowState);
adaptive.elapsed_r3_s=toc(ticFit);adaptive.r3_method="adaptive";

cfgFixed=cfg;
cfgFixed.route30GapHalfWidthMm=0;
cfgFixed.route30GapCount=1;
cfgFixed.route30SingleGapHalfWidthMm=0;

stateFixed=C0.lowState;stateFixed.gHat=C0.pathCal.g0;
ticFit=tic;
fixed=run_single_near_sync_voltage_vp(highMap,C0.templateModel,cfgFixed,stateFixed);
fixed.elapsed_r3_s=toc(ticFit);fixed.r3_method="fixed";

% g_used is the calibrated family's absolute physical-clearance coordinate.
% State-matched therefore fixes the supplied operating state itself. With a
% noisy low calibration this diagnostic is not guaranteed to be a strict
% numerical upper bound, and is reported as a benchmark rather than an oracle.
stateMatched=C0.lowState;stateMatched.gHat=gTruth;
ticFit=tic;
matched=run_single_near_sync_voltage_vp(highMap,C0.templateModel,cfgFixed,stateMatched);
matched.elapsed_r3_s=toc(ticFit);matched.r3_method="state_matched";

Fits=struct('fixed',fixed,'adaptive',adaptive,'state_matched',matched);
end
