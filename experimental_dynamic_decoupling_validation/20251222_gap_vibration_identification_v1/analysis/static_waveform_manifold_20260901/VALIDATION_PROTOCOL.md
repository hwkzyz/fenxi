# Static-to-Dynamic Validation Protocol

This protocol is the regression gate for every change to the gap-library,
low-speed localization, or high-speed migration method. Raw data and formal
legacy outputs are never overwritten; every run writes to an analysis result
directory with a method-specific suffix.

## Fixed test hierarchy

### V1: Same blade, known transition

Use one measured static state as the low-speed anchor and another measured
state of the same blade as the truth after a prescribed gap change. The
method predicts the second waveform from the first. Report waveform RMSE,
peak-voltage error, peak-position error, and gap-change error.

### V2: Same-blade leave-one-gap-out

Remove one nominal gap state from a blade, fit using the remaining states,
and predict the removed waveform. Use the fixed nominal gap labels from the
static experiment. Report per-gap and per-blade error distributions.

### V3: Leave-one-blade-out transfer

Hold out an entire blade. Train the increment/path model on the remaining
blades, use one held-out state as the low-speed anchor, and predict another
held-out state. Absolute waveform mismatch and increment-transfer error must
be reported separately. This is mandatory for any method claiming
cross-blade transfer; it is an external stress test for a method explicitly
restricted to one target blade.

### V4: Experimental low-speed to high-speed closure

Use the actual rotating low-speed template and the corresponding high-speed
windows. Compare reconstructed no-vibration waveforms and final sliding-
window identification against the no-gap baseline and independent vibration
evidence where available.

## Required controls

- Keep the same B2 reference coordinate, OPR alignment, waveform domain,
  valid-support mask, and nominal gap vector for all method variants.
- Never independently peak-align non-reference blades before a physical
  comparison or increment calculation.
- Keep `b_s` out of the formal high-speed increment model. Freeze `a_s` after
  low-speed calibration unless the experiment is explicitly an ablation.
- Separate absolute-template error from gap-increment error.
- Save per-state, per-blade, per-sensor results, not only pooled means.

## Pass/fail record

Each method change must produce:

1. a method identifier and source commit/date;
2. the four validation tables above;
3. per-window high-speed results when V4 is run;
4. the exact data files and configuration used;
5. a short conclusion stating which claims passed, failed, or were not in
   scope.

Pooled averages are not sufficient to declare success. A method is accepted
only for the scope supported by its validation: V1/V2 for same-blade
continuity, V3 for cross-blade transfer, and V4 for the final experimental
claim.

## Current regression-gate thresholds

The supplied gate uses the 90th percentile of per-case errors, so a single
pooled mean cannot hide failed states. The default thresholds are: V1 and V2
relative waveform RMSE <= 10%; and V4 all windows must retain the same EO with
mean absolute frequency difference <= 0.05 Hz, amplitude difference <= 0.001
mm, and RMSE difference <= 0.10 mV. For V3, the primary endpoint is the direct
held-out target waveform error, `RMSE_transfer` (with its P50/P90 reported for
every method). A numerical absolute V3 gate is intentionally not hard-coded
until a repeated-measurement noise/repeatability reference is available. The
ratio `RMSE_transfer/RMSE_zero_increment` and the fraction of cases improved
over zero increment remain secondary diagnostic endpoints. These are explicit
engineering gates for this analysis dataset, not claims of universal
instrument accuracy; any future numerical threshold must be recorded with its
repeatability basis and followed by a rerun of all validation levels.
