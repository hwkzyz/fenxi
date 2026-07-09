# 20251222 BTT Data Foundation

This folder is the first-round foundation scaffold for the 20251222 BTT data
route. Its purpose is to separate reusable data products from method-specific
identification scripts.

## Current Scope

Round 1 only builds a low-risk bootstrap layer:

1. `BTTDataConfig_20251222.m`
   - Stores stable data facts, paths, sensor IDs, OPR channel, blade count,
     sample rate, and output folders.
2. `BTTProjectConfig_20251222.m`
   - Thin project wrapper around the data config.
3. `Step01_Build_LowSpeed_Reference_20251222.m`
   - Adapts the existing 20251222 low-speed `Sensor_Config` product into the
     foundation output layout.
   - Does not re-run raw low-speed extraction yet.

The old route remains untouched:

```text
../20251222_low_speed_rotating_calibration/
```

## First Command

Run from MATLAB inside this folder:

```matlab
Step01_Build_LowSpeed_Reference_20251222
```

Expected outputs:

```text
output/step01_low_speed_reference/
  Sensor_Config_20251222.mat
  Sensor_Config_Summary_20251222.csv
  Standard_Relative_Angles_Quality_20251222.csv
  Step01_Path_Audit_20251222.csv
  LowSpeed_Reference_Metadata_20251222.mat
```

Expected diagnostic figures:

```text
output/figures/step01_low_speed_reference/
  Step01_LowSpeed_Fingerprint_Heatmap_20251222.png
  Step01_LowSpeed_Fingerprint_Heatmap_20251222.pdf
  Step01_Standard_Angle_Quality_20251222.png
  Step01_Standard_Angle_Quality_20251222.pdf
  Step01_Numbering_Overview_20251222.png
  Step01_Numbering_Overview_20251222.pdf
```

These figures are intentionally lightweight in Round 1 because this step is
still a bootstrap adapter. They audit the existing low-speed reference instead
of re-reading raw low-speed waveforms.

## Bootstrap Sources

Step01 checks these sources in order:

```text
../20251222_low_speed_rotating_calibration/output/new_flow/00_low_speed_sensor_config/Sensor_Config_20251222.mat
../20251222_low_speed_gap_prior_decoupling/legacy/output/1000rpm无振动_reference/Sensor_Config_20251222.mat
```

## Step01 Data Contract

`Sensor_Config_20251222.mat` is the main reusable low-speed reference product.
It contains:

```text
Sensor_Config.Dataset
Sensor_Config.ReferenceCase
Sensor_Config.Sensor_IDs
Sensor_Config.Candidate_Sensor_IDs
Sensor_Config.OPR_ID
Sensor_Config.Blades_Num
Sensor_Config.Sample_Rate_Hz
Sensor_Config.Pinlv
Sensor_Config.R_Tip_mm
Sensor_Config.Fingerprints
Sensor_Config.Target_Indices
Sensor_Config.Standard_Relative_Angles
Sensor_Config.Standard_Relative_Angles_OPRCenter
Sensor_Config.Standard_Relative_Angles_Reference
Sensor_Config.OPR_Timing_Method
Sensor_Config.Probe_Arrival_Method
Sensor_Config.Bootstrap_Source_File
Sensor_Config.Bootstrap_Source_Label
```

Additional variables saved in the same file:

```text
SensorConfigSummary
OPRCenterAngleTable
AngleQualityTable
PathAuditTable
Step01Meta
```

`LowSpeed_Reference_Metadata_20251222.mat` is a lighter provenance layer. It
states that Round 1 is `bootstrap_from_existing_sensor_config` and records the
selected source file.

## Path Audit

Step01 writes `Step01_Path_Audit_20251222.csv`, checking:

```text
dataset_root
low_speed_case_dir
dynamic_case_dir
strain_root
newflow_sensor_config
legacy_sensor_config
selected_source
```

The external raw-data folders are audited for reproducibility, but the Round 1
Step01 can still succeed as long as a valid bootstrap `Sensor_Config` exists.

## Difference From 20250527

The route follows the same foundation idea as `20250527_btt_data_foundation`,
but the first round is deliberately narrower:

```text
20250527:
  Step01 re-extracts low-speed BTT features from raw data.

20251222 Round 1:
  Step01 adapts an already accepted NewFlow/legacy Sensor_Config.
```

Current 20251222 experiment facts:

```text
sensor_ids = [1 2 3]
candidate_sensor_ids = [1 2 3]
opr_id = 4
blades_num = 6
sample_rate_hz = 5e6
r_tip_mm = 62.0
low_speed_case = 1000rpm无振动
dynamic_cases = {1000_2500_3500}
```

Future rounds should replace the bootstrap adapter with raw low-speed feature
extraction only after Step02/Step03 interfaces are stable.

## Step02 Dynamic Raw Extraction

Round 2 is now implemented by `Step02_Extract_Dynamic_BTT_20251222.m`.

Run from MATLAB inside this folder:

```matlab
Step02_Extract_Dynamic_BTT_20251222
```

Step02 reads the raw dynamic MAT files directly from:

```text
cfg.dataset_root / 1000_2500_3500 / 4-<channel>-<file_id>.mat
```

It does not read old dynamic `jiluOPR.mat`, `omega.mat`, or
`jilublade_probe*.mat` files from NewFlow or legacy folders.

Produce standard method-neutral dynamic products:

```text
output/step02_dynamic_btt/1000_2500_3500/
  DynamicBTTFeature_20251222.mat
  DynamicPulseTable_20251222.csv
  SensorNumberingQuality_20251222.csv
  RevolutionNumberingQuality_20251222.csv
  OPRTable_20251222.csv
  OmegaTable_20251222.csv
  OPRSemanticsCheck_20251222.csv
  Step02_SourceInventory_20251222.csv
  Step01_Step02_Consistency_20251222.csv
  Step02_Summary_20251222.csv
  Step02_Metadata_20251222.mat
```

Also export newly generated compatibility files for comparison and temporary
reuse inside this foundation output folder:

```text
output/step02_dynamic_btt/1000_2500_3500/
  jiluOPR.mat
  omega.mat
  jilublade_probe1.mat
  jilublade_probe2.mat
  jilublade_probe3.mat
```

Step02 is a data layer only. It does not run direct-template identification,
gap-prior decoupling, eta/dx fitting, VP screening, or TARC validation.

Important Step02 rules:

```text
jiluOPR(:,1) = opr_center_time_s
jiluOPR(:,2) = opr_start_time_s
jiluOPR(:,3) = opr_end_time_s
OPR centers are computed from the raw OPR waveform by multi-threshold width centers
opr_events_per_revolution is inferred from raw OPR centers and expected RPM range
rev_id is assigned by arrival time relative to OPRTable
blade_id is assigned from same-revolution feature vectors matched to Step01 fingerprints
assignment_quality is the same-revolution fingerprint correlation, not a validity flag
```

For the current `1000_2500_3500` raw case, Step02 selected:

```text
opr_events_per_revolution = 6
OPR raw pulse count = 37845
dynamic pulse count = 113535
valid dynamic pulse count = 113526
rpm_median = about 2922.9
sensor quality = good for sensors 1, 2, and 3
```

`Step02_SourceInventory_20251222.csv` records the raw channel file inventory:

```text
channel 1, 2, 3 = BTT probe channels
channel 4       = OPR channel
```

Expected diagnostic figures:

```text
output/figures/step02_dynamic_btt/1000_2500_3500/
  Step02_RPM_Trend_20251222.png
  Step02_OPR_Center_Check_20251222.png
  Step02_BladeID_Sequence_20251222.png
  Step02_Revolution_Coverage_20251222.png
  Step02_Fingerprint_Match_20251222.png
  Step02_Pulse_Feature_Trend_20251222.png
```

## Step03 Dynamic Observation Bundle

Round 3 is now implemented by `Step03_Build_Dynamic_Observation_Bundle_20251222.m`.

Run from MATLAB inside this folder:

```matlab
Step03_Build_Dynamic_Observation_Bundle_20251222
```

Create the reusable observation layer:

```text
output/step03_observation_bundle/1000_2500_3500/
  BTT_Observation_Bundle_20251222.mat
  BTT_Observation_LongTable_20251222.csv
  BTT_Observation_Summary_20251222.csv
  BTT_Observation_QualityFlags_20251222.csv
  Step03_AngleSourceCheck_20251222.csv
  Step03_RevIDConsistency_20251222.csv
  Step03_OmegaConsistency_20251222.csv
  Step03_InputConsistency_20251222.csv
  Step03_Metadata_20251222.mat
  jilublade_probe1_vib_final.mat
  jilublade_probe2_vib_final.mat
  jilublade_probe3_vib_final.mat
```

Step03 consumes only the new foundation products:

```text
Step01: output/step01_low_speed_reference/Sensor_Config_20251222.mat
Step02: output/step02_dynamic_btt/1000_2500_3500/DynamicBTTFeature_20251222.mat
```

It converts each dynamic pulse row into:

```text
observed_angle_deg      = Step02 dynamic relative angle
standard_angle_deg      = Step01 low-speed OPR-center standard angle
angle_deviation_deg     = wrapped(observed_angle_deg - standard_angle_deg)
displacement_mm         = angle_deviation_deg * pi/180 * r_tip_mm
rpm                     = recomputed from Step02 OPRTable
```

Step03 also verifies the Step02 semantics before downstream use:

```text
Step03_InputConsistency_20251222.csv     checks Step01/Step02 sensor, blade, OPR, and timing metadata
Step03_AngleSourceCheck_20251222.csv     compares Step02 relative_angle_deg with OPRTable fallback angle
Step03_RevIDConsistency_20251222.csv     compares Step02 rev_id with arrival_time + OPRTable rev_id
Step03_OmegaConsistency_20251222.csv     compares OPRTable-recomputed RPM with Step02 OmegaTable
```

For the current `1000_2500_3500` case, Step03 produced:

```text
observation rows = 113535
valid rows       = 113517
valid ratio      = about 99.984%
relative vs fallback angle median difference = 0 deg
rev_id mismatch ratio = 0
OPRTable RPM vs OmegaTable RPM median difference = 0 rpm
```

Expected diagnostic figures:

```text
output/figures/step03_observation_bundle/1000_2500_3500/
  Step03_Displacement_Trend_20251222.png
  Step03_AngleDeviation_Trend_20251222.png
  Step03_RPM_And_Coverage_20251222.png
  Step03_Displacement_Distribution_20251222.png
  Step03_Quality_Flags_20251222.png
  Step03_Observed_vs_Standard_Angle_20251222.png
  Step03_AngleSource_Check_20251222.png
  Step03_RevID_Mismatch_20251222.png
```

Step03 is still a data layer. It does not run template identification,
gap-prior decoupling, eta/dx fitting, VP screening, or TARC validation.

## Step04 Low-Speed Template Library

Round 4 is now implemented by `Step04_Build_LowSpeed_OPRCenterStd_Template_20251222.m`.

Run from MATLAB inside this folder:

```matlab
Step04_Build_LowSpeed_OPRCenterStd_Template_20251222
```

Step04 consumes the Step01 `Sensor_Config_20251222.mat`, but it does not use
old low-speed template files from the legacy route. It directly re-reads the
raw low-speed case:

```text
cfg.dataset_root / cfg.low_speed_case / 4-<channel>-<file_id>.mat
```

It extracts low-speed OPR centers with the same multi-threshold center method
used by Step02, extracts CH1/CH2/CH3 probe pulse windows, maps raw waveform
samples to the OPRCenterStd spatial coordinate, and builds all 18
sensor-blade templates:

```text
3 sensors x 6 blades = 18 templates
```

Create the accepted low-speed template library products:

```text
output/step04_low_speed_template/
  Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat
  Template_OPRCenterStd_Summary_20251222.csv
  LowSpeed_Template_SourceData_20251222.mat
```

Expected diagnostic figures:

```text
output/figures/step04_low_speed_template/
  Step04_Template_Curves_Clean_CH1.png
  Step04_Template_Curves_Clean_CH2.png
  Step04_Template_Curves_Clean_CH3.png
  Step04_RawPointCloud_TemplateFit_CH1.png
  Step04_RawPointCloud_TemplateFit_CH2.png
  Step04_RawPointCloud_TemplateFit_CH3.png
  Step04_Template_Quality_Heatmap.png
  Step04_Template_Summary_Table.png
```

For the current low-speed case, Step04 produced:

```text
low-speed OPR pulses = 837
CH1 low-speed probe pulses = 836
CH2 low-speed probe pulses = 836
CH3 low-speed probe pulses = 837
template entries = 18
good/usable template entries = 18
```

## Step05: Full-Raw Direct-Template Identification

The only official Step05 entry point is
`Step05_SingleSync_DirectTemplate_Identification_20251222.m`.

Run from MATLAB inside this folder:

```matlab
Step05_SingleSync_DirectTemplate_Identification_20251222
```

This script follows the `20250527_btt_data_foundation` Step05 route. It does
not read products from `20251222_low_speed_rotating_calibration`, and it does
not use `pulse_feature_value -> template displacement` as the identification
observation. Its main data chain is:

```text
Step02 DynamicBTTFeature_20251222.mat
  DynamicPulseTable supplies source_file_id, pulse windows, OPRTable

Step03 BTT_Observation_Bundle_20251222.mat
  Observation_Table supplies blade/rev/sensor labels and valid observations

Step04 Template_OPRCenterStd_LowSpeed_AllBlades_S123_20251222.mat
  low-speed template waveform T_{sensor,blade}(x)

Raw dynamic waveform files
  cfg.dataset_root / 1000_2500_3500 / 4-<sensor>-<file_id>.mat
```

Forward model:

```text
V_i ~= T_{sensor,blade}(x_i - dx_c - eta_s - A*sin(EO*theta_i + phi))
```

The method is the 20250527 full-raw query-safe route adapted for the
unconstrained 20251222 frequency-audit requirement:

```text
1. Re-read raw dynamic waveform windows around blade passages.
2. Map every raw sample into OPRCenterStd x coordinates.
3. Build a query-safe waveform bundle against the Step04 template domain.
4. Evaluate all finite EO candidates in the search band; no EO prior and no
   VP top-k pruning are used by default.
5. Refine A, phi, and dx_c by full waveform shape-correlation objective.
   Default sensor eta fitting is disabled.
6. Use multi-window joint EO evidence for the final frequency decision.
   Ambiguous joint evidence is kept as diagnostics and does not produce a
   final frequency.
7. Save trend, reconstruction previews, masks, EO landscape, and joint EO
   diagnostics.
```

Current official default settings are intentionally unconstrained:

```text
objective_mode = corr
top_k_eo = Inf
amplitude_limit_mm = Inf
dx_c_limit_mm = Inf
sensor_eta_limit_mm = Inf
sensor_eta_reg_weight_v_per_mm = 0
overshoot_penalty_weight = 0
joint_eo_gap_ratio_threshold = 0.03
```

`invalid_query_penalty_scale` is only a template-domain validity protection:
it prevents a candidate from gaining score by moving query points outside the
calibrated Step04 template domain. It is not an EO prior and not a physical
parameter bound.

Main outputs:

```text
output/step05_single_sync_direct_template/1000_2500_3500/
  Result_Step05_FoundationFullRawTop3RMSE_B<blade>_S123_20251222.mat
  Trend_Step05_FoundationFullRawTop3RMSE_B<blade>_S123_20251222.csv

output/step05_single_sync_direct_template/
  Step05_TargetBlade_Summary_20251222.csv
  Step05_RunInfo_20251222.mat
```

Useful test-run environment variables:

```text
STEP05_TARGET_BLADES=1
STEP05_TARGET_LAPS=3
STEP05_WINDOW_LAPS=3
STEP05_SLIDING_STEP_LAPS=1
STEP05_DEBUG_MAX_WINDOWS=1
STEP05_OUTPUT_TAG=smoke_w3
STEP05_SHOW_PLOTS=false
STEP05_SAVE_FIGURES=false
```

Recommended Step05 validation order:

```text
1. Single target blade + 2 windows:
   STEP05_TARGET_BLADES=1
   STEP05_DEBUG_MAX_WINDOWS=2

2. Same target blade + all target windows:
   STEP05_TARGET_BLADES=1
   STEP05_DEBUG_MAX_WINDOWS=inf

3. Run additional comparison blades one by one only when needed:
   STEP05_TARGET_BLADES=4
```

Avoid defaulting to a full-blade sweep while the Step05 method is being
cleaned up. Blade-specific template quality, query-safe coverage, and dynamic
waveform quality can otherwise obscure method-logic problems.

Keep `STEP05_OUTPUT_TAG` short, preferably no more than 24 characters. Long
tags can make MATLAB `-v7.3` result-file saves fail on Windows path-length
limits.

Current frequency audit evidence:

```text
output/diagnostics/Step05_Unconstrained_FrequencyEvidence_B1B2B3_20251222.csv
```

Run the policy audit from MATLAB inside this folder:

```matlab
Audit_Step05_UnconstrainedFrequencyPolicy_20251222
```

Expected audit outputs:

```text
output/diagnostics/Step05_UnconstrainedPolicy_Audit_20251222.csv
output/diagnostics/Step05_UnconstrainedFrequency_ResultAudit_20251222.csv
output/diagnostics/Step05_UnconstrainedFrequency_OverallAudit_20251222.csv
```

Summary of the audited blades:

```text
B1: legacy EO14, foundation joint EO14, final frequency accepted.
B2: legacy EO7, foundation joint evidence ambiguous, final frequency withheld.
B3: legacy EO24, foundation joint evidence ambiguous, final frequency withheld.
```

This is intentional: the official Step05 should output a frequency only when
the unconstrained waveform evidence is strong enough. Otherwise it leaves the
final frequency as `NaN` and preserves the EO candidate table for diagnosis.

## Optional Diagnostics

The feature-level template-guided baseline has been moved out of the main
Step sequence:

```text
optional_diagnostics/
  Diag01_TemplateGuided_FeatureBaseline_20251222.m
```

It is not the official Step05, is not used by default downstream steps, and
does not produce the full-raw direct-template result. It only uses
`pulse_feature_value` plus the Step04 template to build a feature-level
baseline for comparison.

Run only when a diagnostic comparison is needed:

```matlab
run(fullfile('optional_diagnostics', 'Diag01_TemplateGuided_FeatureBaseline_20251222.m'))
```

Diagnostic outputs are intentionally separated from the Step05 directory:

```text
output/optional_diagnostics/template_guided_feature_baseline/1000_2500_3500/
  Diag01_TemplateGuided_FeatureBaseline_Result_20251222.mat
  Diag01_TemplateInversion_Observation_20251222.csv
  Diag01_TemplateGuided_FeatureBaseline_Trend_20251222.csv
  Diag01_TemplateGuided_FeatureBaseline_Summary_20251222.csv
  Diag01_Metadata_20251222.mat

output/figures/optional_diagnostics/template_guided_feature_baseline/1000_2500_3500/
  Diag01_Amplitude_Trend_20251222.png
  Diag01_EO_RMSE_20251222.png
  Diag01_TemplateInversion_Check_20251222.png
```

After these data layers are stable, method-specific identification scripts
such as direct-template with `eta_s`, gap-prior routes, or VARPRO variants can
consume the standard outputs instead of rebuilding their own data inputs.
