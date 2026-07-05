# 20241106 BTT Data Foundation

This folder refactors the existing `20241106_low_speed_rotating_calibration`
route into the same data-foundation style used by `20250527_btt_data_foundation`.

## Target Workflow

```matlab
Step01_Build_LowSpeed_Reference_20241106
Step02_Extract_Dynamic_BTT_20241106
Step03_Build_Dynamic_Observation_Bundle_20241106
Step04_Build_LowSpeed_OPRCenterStd_Template_20241106
Step05_SingleSync_DirectTemplate_Identification_20241106
Step06A_Build_StrainRPM_ResonanceEvidence_20241106
Step06B_Validate_Step05_With_Strain_TARC_20241106
```

## Round 1 Scope

Round 1 creates the data-foundation scaffold:

- `BTTDataConfig_20241106.m`
- `BTTProjectConfig_20241106.m`
- this `README.md`

The config intentionally stores only data facts and extraction/reference
settings. Direct-template, gap-library, eta, dx, EO-search, and optimizer
settings should stay in later method scripts.

## Round 2 Scope

Round 2 builds the low-speed reference layer from the original low-speed
rotating-calibration logic. It writes:

```text
output/step01_low_speed_reference/
  Sensor_Config_20241106.mat
  LowSpeed_Features_20241106.mat
  LowSpeedReference_20241106.mat
  LowSpeedReferenceFingerprint_20241106.mat
  LowSpeedNumbering_20241106.mat
  Sensor_Config_Summary_20241106.csv
  LowSpeed_Channel_Stats_20241106.csv
  LowSpeedNumbering_20241106.csv
  LowSpeedNumbering_Summary_20241106.csv
  Standard_Relative_Angles_Quality_20241106.csv
```

## Design Rule

Compatibility files may be emitted for old scripts, but old outputs must not
be treated as the source of truth for new data-foundation steps. Step02 should
independently extract dynamic BTT timing from the raw `3000_3150` case.

Step01 may use start-edge OPR times internally for low-speed revolution
segmentation, but the exported standard coordinate is OPR-center based:
`Standard_Relative_Angles_OPRCenter` is the downstream angle source of truth.
Step04 reads `Sensor_Config.OPR_Center_Times` for low-speed template
coordinates.

## Round 3 Scope

Round 3 extracts dynamic BTT timing from the raw `3000_3150` case. The new
source-of-truth product is:

```text
output/step02_dynamic_btt/3000_3150/
  Step02_Dynamic_BTT_Extraction_20241106.mat
  Step02_Summary_20241106.csv
```

Legacy-compatible files are also emitted for adapter scripts:

```text
  jiluOPR.mat
  omega.mat
  jilublade_probe2.mat
  jilublade_probe3.mat
  jilublade_probe5.mat
  jilublade_probe7.mat
```

Step02 exports canonical dynamic OPR timing with `multi_threshold_center` in
`jiluOPR(:,1)` and keeps start/end edge columns only for audit compatibility.
Probe arrival timing uses `polynomial_peak_fit`. The dynamic OPR center axis
has been checked against the original NewFlow/legacy OPR file; the median
difference is 0 s and the maximum absolute difference is about `8.1e-08` s.
On the full run, all four analysis sensors passed dynamic blade labeling with
`quality_status=good`.

## Step03 Observation Bundle

Step03 builds the method-neutral arrival-time observation layer. It uses
`opr_pulses_per_rev = 1`, `Sensor_Config.Standard_Relative_Angles_OPRCenter`,
and the Step02 dynamic labels. It must not store raw waveform slices or run
template fitting; method steps read raw data themselves.

## Step04-Step06 Scope

The later method steps follow the `20250527_btt_data_foundation` structure:
configuration-driven inputs, source checks, metadata, summary tables, quality
flags, restartable outputs, and diagnostic figures.

The physical definitions follow `20241106_low_speed_rotating_calibration`:
analysis sensors `[2 3 5 7]`, one OPR pulse per revolution, low-speed
OPR-center standard angles, and direct-template identification that reads
dynamic raw waveform slices at the method layer.

## Step04 Two-Stage Calibration

Step04 uses a two-stage low-speed waveform calibration policy:

```text
wide pulse-window point cloud
  -> preliminary wide template
  -> trusted x_domain selection
  -> raw points inside x_domain only
  -> re-bin / median / smooth
  -> final x_grid / v_grid calibration template
```

The complete wide point cloud is retained for waveform-shape inspection,
template-center checks, coverage checks, trusted-domain selection, and
diagnostic figures. It is not used wholesale as the final calibration
template.

Each `Template.SensorBlade` entry now separates diagnostic and formal
calibration fields:

```text
x_wide_domain              full collected pulse-window x range
x_domain                   trusted calibration domain
x_query_safe_domain        x_domain shrunk by the Step04 query-safe margin
wide_x_grid / wide_v_grid  preliminary wide-domain diagnostic template
x_grid / v_grid            final template re-aggregated from trusted points
final_template_point_policy = trusted_domain_recalibration
```

The saved `point_cloud` entries include `is_inside_x_domain` and
`is_used_for_final_template` so figures can distinguish wide waveform points
from the low-speed points that actually generated the final template.

## Target Blade Rule

Do not infer the target blade from convenience defaults. The active 20241106
rotating-calibration identification script sets:

```matlab
S06.blades = 4;
S06.outputLabel = 'B4_only';
```

in `20241106_low_speed_rotating_calibration/Step06_RunIdentificationByBlade_20241106.m`.
Therefore the default foundation Step05/Step06 validation target is B4.
All-blade results are allowed only as an explicit scan/audit path
(`cfg.step05_scan_blades = 1:cfg.blades_num`), not as the main validation
claim.

## 2026-07-04 Validation Snapshot

The current foundation route was rerun with unified OPR-center timing:

- Step02: `multi_threshold_center`, OPR count `8514`; labeled pulses CH2/CH3/CH5/CH7 = `51084/51078/51084/51078`.
- Step03: observation rows `204350`, valid rows `204324`; RPM is recomputed from `jiluOPR(:,1)` with `opr_pulses_per_rev = 1`.
- Step04: 24/24 low-speed OPRCenterStd templates are good/usable.
- Step05: main target B4 with S257, 17/17 windows `ok`, EO = 12, frequency range `629.143-629.166 Hz`.
- Step06A: one strain/RPM resonance region, dominant order 12, peak frequency about `631.193 Hz`.
- Step06B: 17/17 windows validated, median absolute frequency error about `0.0071 Hz`, mean best TARC about `0.970`.

The legacy/NewFlow OPR file is now used only as an audit reference. The
foundation source of truth for Step05 timing is Step02 `jiluOPR(:,1)`.
