# Route 30: low/high-speed gap and vibration identification

The frozen synchronous dual-order simulation method is documented in
`DUAL_SYNC_FUNNEL_V2_SPECIFICATION.md`. Its formal entry is
`Run_DualSyncFunnelV2Generalization`.

The final-candidate competition audit is documented in
`RESULTS_DUAL_SYNC_FINAL_COMPETITION_20260810.md`. It defines the validated
scope of blind EO selection and the close-order boundary.

The 2026-08-09 dual-frequency efficiency audit is in
`RESULTS_FAST_DUAL_CANDIDATE_20260809.md`. Route 30 now uses a
derivative-weighted conditional batch Top-K candidate generator with an
automatic full two-dimensional fallback for low-confidence support.

The synchronous-order fast-route benchmark is in
`RESULTS_SYNCHRONOUS_ORDER_FAST_20260809.md`. The main entry accepts
`cfg.route30FrequencyStructureMode` so synchronous records can search integer
EO directly instead of entering the general 5 Hz two-frequency grid.

This route is the formal shared-gap simulation main through
`run_inv_log_2_main_method`. The low-speed-template increment model is
available through `cfg.route30ForwardModel="low_increment"` and is kept as a
separate validation mode until its cached VP implementation is fully profiled.
The original scripts in this directory remain validation and regression drivers.

Formal entries:

```matlab
R = run_inv_log_2_low_high_main(lowSpeedInput, highMap, templateLib, cfg);
R = run_inv_log_2_main_method(highMap, templateLib, cfg, lowSpeedState);
```

Optional per-sensor extension:

```matlab
R = run_per_sensor_gap_extension(highMap, templateLib, cfg, lowSpeedState);
```

The extension fits one absolute gap per sensor after the shared frequency
seed. For difficult per-sensor cases, `cfg.perSensorFrequencySeed` may be
provided explicitly; this is an extension diagnostic, not a default prior.

The first form is preferred because it makes the low-speed calibration input
explicit. The second form is retained for compatibility with existing callers.

For `cfg.route30ForwardModel="low_increment"`, a raw low-speed waveform now
uses the promoted low-speed template calibration by default: 0.02 mm spatial
bin means followed by complete-revolution grouped cross-validation over
physical Savitzky--Golay window lengths. The selected window, span, and CV
curve are returned in `R.low_speed_selected_window_mm`,
`R.low_speed_selected_span_bins`, and `R.low_speed_template_cv`. Selection uses
only the low-speed record. Set `cfg.route30LowTemplateMethod="legacy"` only to
reproduce the former per-turn-interpolation aggregation.

For repeated high-speed records under one measurement setup, calibrate once
and reuse the result:

```matlab
C = calibrate_inv_log_2_low_speed(lowSpeedInput, templateLib, cfg);
R1 = run_inv_log_2_high_with_calibration(C, highMap1, cfg);
R2 = run_inv_log_2_high_with_calibration(C, highMap2, cfg);
```

The calibration object stores the selected physical window, CV curve,
low-speed path calibration, sampled forward template, and cached increment
model; it does not retain the raw mapped low-speed samples.

## Method

1. Generate a noisy low-speed, no-vibration waveform from one curve in the gap library.
2. Aggregate turns in the spatial domain and estimate `g_low` and `dx_low` against the full library.
3. Use `g_low_hat` only as the center of a broad high-speed gap profile; `g_high` remains free.
4. Generate the high-speed waveform at `g_high`, using the same SNR definition as the low-speed waveform.
5. Run an independent single-frequency VP scan and a baseline unordered dual-frequency VP scan, then select the provisional model order by BIC.
6. For provisional dual records only, preserve the baseline candidate pool and add one-dimensional conditional scans that anchor one frequency while scanning the other. This avoids a second full two-dimensional VP scan.
7. Reproject every retained dual candidate on the complete waveform, rather than reusing amplitudes and phases estimated on the screening subset.
8. Short-refine the best 10 full-waveform candidates and spend the full nonlinear iteration budget only on the leading basin.
9. Apply query-safe masking only to the finite-support low-speed increment model. The absolute template keeps the complete trusted sample set.

The pointwise inverse map is not used by the main route.

## Entry points

```matlab
Run_LowSpeedGapAudit('smoke')
Run_LowHighVPBridge('smoke')
Run_LowHighVPBridge('full')
Run_LowHighSingleVPBridge('smoke')
Run_LowHighSingleVPBridge('full')
Run_SequentialDualCandidateAudit('smoke')
Run_FormalMainSmoke()
Run_FormalMainRegression()
Run_FormalDualRepeatRegression(1:3)
```

## Verified results

- Unified formal regression (2026-08-07): single-frequency 3/3 and dual-frequency 3/3 successes. The maximum frequency error was 0.092 Hz for the three dual cases and 0.065 Hz for the three single cases.
- Independent dual repeat (2026-08-07): 9/9 successes for `500+1300`, `700+1200`, and `900+1400 Hz` over three high-speed noise seeds and unequal low/high gap directions. Maximum frequency error was 0.201 Hz. High-speed identification took about 11.4--21.4 s per case in this diagnostic run.
- The previously failing `g_low=0.8 mm -> g_high=0.2 mm`, `700+1200 Hz` case had the true pair at full-waveform replay rank 4. Conditional completion plus Top-10 short refinement recovered it for all three noise seeds.

- Dual-frequency full: 30/30 parameter successes; maximum 20 dB frequency error 0.270 Hz.
- Non-grid dual smoke with a 5 Hz VP grid: 6/6 successes; maximum frequency error 0.165 Hz.
- Non-grid dual repeat (3 independent low/high noise seeds, two gap pairs): 18/18 successes; maximum frequency error 0.214 Hz; mean high-speed fitting time about 6.9 s per case.
- Non-grid dual low-SNR audit: 12/12 successes at 15 and 10 dB with 900 samples. At 5 dB, 900 samples gave 4/6, while increasing to 1500 stratified samples restored 6/6 with maximum frequency error 0.823 Hz and about 7.4 s high-speed fitting per case.
- Adaptive 5 dB repeat: start with 1500 samples and rerun with 2500 only when `rhoJ < 0.003`. Across 18 cases and three independent low/high noise seeds, 4 cases used the dense fallback and all 18 succeeded; maximum frequency error was 1.14 Hz and mean total high-speed fitting time was about 13.4 s per case.
- Broad random-frequency audit at 20 dB: six deterministic off-grid pairs across two gap pairs all succeeded (12/12), with maximum frequency error 0.411 Hz and mean high-speed fitting time about 11.9 s per case.
- Random amplitude/phase audit at 20 dB: four off-grid pairs, four amplitude pairs, four phase pairs, and two gap pairs succeeded 8/8 after weak-component fallback; the normal 900-sample pass used the 1500-sample rerun twice.
- Robust dual full with a 5 Hz grid across the original five gap pairs and no/20 dB conditions: 30/30 successes; maximum frequency error 0.275 Hz; about 8.3 s high-speed fitting per case.
- Single-frequency full: 90/90 parameter and BIC order successes; maximum 20 dB frequency error 0.202 Hz.
- Low/high gap pairs include equal and unequal gaps.
- Low- and high-speed records use the same SNR per case.
- Single `mu=0` final refinement reduced dual high-speed fit time from about 100.7 s to 83.9 s and single full wall time from about 333 s to 276 s without changing success counts.

## Rejected shortcuts

- Reducing the high-speed gap grid from 61 to 31 points produced only 6/9 dual smoke successes.
- Sequential first-frequency/residual-second-frequency candidates recalled only 3/6 dual smoke truths.
- Per-gap local frequency refinement increased runtime.
- A 12.5 Hz non-grid grid achieved only 5/6 successes because a narrow true-frequency basin was crowded out; the 5 Hz grid restored 6/6 at about 7.2 s of high-speed fitting per case.

These alternatives remain diagnostic only and are not enabled in the verified route.

## Structured four-mode entry

The formal known-structure entry is now:

```matlab
R = run_inv_log_2_structured_main(lowSpeedInput, highMap, templateLib, cfg, mode);
```

Supported modes are `single_sync`, `single_async`, `dual_sync_sync`, and
`dual_sync_async`. Synchronous components use the measured shaft-angle basis;
asynchronous components use continuous time-frequency bases. Every trusted
sample is retained. The default information requirement is eight observed
revolutions, and shorter records are rejected with
`invlog2:InsufficientRevolutions` rather than silently forced to a unique answer.

Each result reports both the selected solution and its nearest distinct refined
competitor through `frequency_margin`, `noise_normalized_margin`,
`identification_confident`, and `competing_frequency_or_order`. Parameter
correctness in a simulation and statistical separation of the winning basin are
therefore reported separately.

Unified regression:

```matlab
Run_StructuredFourModeRegression
Run_StructuredFourModeRepeatRegression([0 1000 2000])
```

See `RESULTS_STRUCTURED_FOUR_MODE_20260807.md` for the 5 dB, eight-revolution
acceptance run and its interpretation.

The dual-mode defaults retain a wider candidate pool (`100` candidates per
gap, followed by up to `2100` iterated voltage replays) than the original
diagnostic setting. This is required for weak components near `0.1 mm`; a
small VP Top-K can remove their basin before the nonlinear voltage stage.
The wider pool costs time but does not remove trusted samples or alter the
final voltage objective.

`structuredMinRevolutions=8` is an execution minimum, not a uniqueness proof.
When the refined objective has a distinct competitor without a sufficient
noise-normalized margin, the result is marked `identification_status="ambiguous"`
and `requires_more_observations=true` (default confidence target: 16 turns).

Boundary evidence and the selective-reporting rule are documented in
`RESULTS_STRUCTURED_BOUNDARY_20260807.md`.

For the declared `A >= 0.1 mm` scope, dual-mode final selection uses a
noise-tolerant component floor of `0.075 mm`. This removes numerical
two-component fits whose second component is physically below the target
scope; it is not used to force a frequency when two admissible basins remain.

The weak-component sensor-layout audit is documented in
`RESULTS_STRUCTURED_SENSOR_LAYOUT_20260807.md`. It shows that candidate and
optimizer changes alone cannot remove all 5 dB aliases from the current
three-sensor layout; a five-sensor nonuniform layout passed the complete
12-case boundary matrix and the two hardest cases over three noise seeds.

Verified settings can be applied explicitly without changing the measured
hardware entry:

```matlab
cfg = configure_structured_identification(cfg,"experiment_current");
cfg = configure_structured_identification(cfg,"simulation_robust5");
```

The second profile is an expanded simulation/design configuration. It must not
be used to reinterpret data collected by the existing three physical probes.

## Paper-scope result

The paper-facing validation uses only the existing three-probe layout, eight
revolutions, all trusted samples, amplitudes at or above `0.1 mm`, and
`10/15/20 dB`. Across the four structured modes and two additional weak-dual
cases, all 42 tests passed; the maximum frequency error was `0.0935 Hz` and the
maximum absolute gap error was about `0.00062 mm`. See
`RESULTS_STRUCTURED_PAPER_SCOPE_20260807.md`.

The 5 dB boundary and five-probe studies are extensions only and are not part
of the paper's main correctness claim.
