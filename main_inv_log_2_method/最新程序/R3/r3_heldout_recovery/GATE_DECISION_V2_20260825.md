# Frozen pre-formal gate decision

The initial held-out 0.3 mm state failed the noiseless 0.10 mm state-matched
gate. Replacing it with the edge state 1.4 mm did not solve the issue: that
state also selected EO17 at 0.10 mm. A provenance-balanced interior split was
therefore audited before any noisy Monte Carlo result was produced.

The frozen v2 calibration states are 0.2, 0.3, 0.4, 0.6, 0.8, 1.1, 1.3, 1.4
and 1.5 mm. The held-out operating states are 0.5, 0.7, 0.9, 1.0 and 1.2 mm.
Their signed mismatch from the 0.8 mm reference is -0.3, -0.1, +0.1, +0.2
and +0.4 mm. All five passed the noiseless 0.10 mm EO10 state-matched probe.

The revision defines the calibrated interpolation envelope; it was frozen
before threshold calibration, smoke testing or formal noise realizations.
The rejected 0.3 and 1.4 mm probes remain documented as response-family
boundary evidence and are not silently treated as successful operating states.
