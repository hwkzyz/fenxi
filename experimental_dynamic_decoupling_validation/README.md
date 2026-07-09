# Experimental Dynamic Decoupling Validation

## Goal

Use the two existing rotor experiments to check whether the proposed
clearance-vibration decoupling method can be applied to real dynamic
waveforms before planning new hardware tests.

This folder is independent from the previous super-Gaussian experiment
workflow. The old workflow is used only as a data extraction and waveform
windowing source; the identification stage should be replaced by the
current joint clearance-vibration model.

## Available Data

### 20251222 experiment

- Raw BTT data:
  `E:\试验数据\20251222\传感器数据`
- Strain-gauge data:
  `E:\试验数据\20251222\应变片数据`
- Existing processing scripts:
  `E:\0小论文+程序\0博士期间小论文+程序\7超高斯模型-权重-瞬态\程序\直叶片验证\实验验证\超高斯\20251222适配-3号传感间隙变化-改`
- Main cases:
  - Low-speed reference: `1000rpm无振动`
  - Dynamic case: `1000_2500_3500`
- Sensors:
  - Capacitance probes: `1, 2, 3`
  - OPR: `4`

This dataset is the preferred first validation target because the sensor
set is cleaner and the dynamic case is long enough for windowed analysis.

### 20250527 experiment

- Raw BTT data:
  `E:\试验数据\20250527\试验20250527`
- Strain-gauge data:
  `E:\试验数据\20250527\应变片数据20250527`
- Existing processing scripts:
  `E:\0小论文+程序\0博士期间小论文+程序\7超高斯模型-权重-瞬态\程序\直叶片验证\实验验证\超高斯\20250527适配`
- Main cases:
  - Low-speed reference: `20250526_910`
  - Dynamic cases: `20250526_3150`, `20250526_2500-3500_t400`
- Sensors:
  - Capacitance probes: `1, 3, 6`
  - Eddy-current probes: `7, 8`
  - OPR: `4`

This dataset is better suited as a second validation target after the
method runs on 20251222. The first pass should use only capacitance probes
`1, 3, 6`; eddy-current probes can be kept for a separate robustness check.

## Core Idea

The old super-Gaussian workflow already solves the difficult experimental
data problem:

1. Read exported `4-channel-block.mat` files.
2. Extract OPR pulses and rotor speed.
3. Segment probe pulses.
4. Match blade IDs using low-speed fingerprints.
5. Extract the same blade across multiple revolutions.
6. Map each pulse window into a local spatial coordinate.

The proposed decoupling method needs a waveform-domain input rather than a
single arrival time. Therefore, the key adaptation is:

```matlab
experiment_case.Extracted_Data(sensor).Laps
    -> highMap.t_v
    -> highMap.x_v
    -> highMap.V_a
    -> highMap.rev_v
    -> highMap.S_v
```

After this conversion, the existing method entry can be called:

```matlab
result = run_inv_log_2_main_method(highMap, templateLib, cfg, staticState);
```

## Validation Questions

The first-pass analysis should answer four practical questions.

1. Can real dynamic pulse windows be converted into the same `highMap`
   structure used by the proposed method?
2. Does the joint model reduce waveform residuals compared with a fixed
   clearance template?
3. Are the identified vibration frequencies consistent with the independent
   strain-gauge spectrum?
4. Does estimated clearance remain physically reasonable across sensors,
   blades, and time windows?

At this stage, the strain-gauge data should be treated as an independent
frequency reference, not as an absolute vibration displacement truth.

## Proposed Workflow

### Step 0: Isolate reusable paths

Create dataset-specific config files in this folder:

- `build_expdec_config_20251222.m`
- `build_expdec_config_20250527.m`

Each config should store:

- raw data path;
- old extraction script path;
- output path;
- selected dynamic case;
- selected sensors;
- target blades;
- start time and lap count;
- template-library settings;
- strain-gauge reference settings.

### Step 1: Reuse old Step1 and Step2

For each dataset, run or load the old outputs:

- `Sensor_Config_*.mat`
- `jiluOPR.mat`
- `omega.mat`
- `jilublade_probe*.mat`

These outputs provide blade identity, OPR timing, and target pulse windows.

### Step 2: Build compact experiment cases

Reuse the old functions:

- `prepare_identification_case_single_sync_20251222.m`
- `prepare_identification_case_single_sync_20250527.m`

The output `experiment_case` contains the real waveform samples for each
selected blade passage:

- `t_points`: sample time;
- `x_points`: local spatial coordinate;
- `v_points`: measured voltage;
- `lap_id`: revolution/window index.

### Step 3: Convert experiment case to highMap

Write a new adapter:

- `experiment_case_to_highmap.m`

The adapter should:

- concatenate selected sensors and laps;
- keep only the trusted local spatial window;
- remove obvious baseline or isolated invalid samples only when necessary;
- assign `rev_v` from lap/window ID;
- assign `S_v` from sensor ID;
- preserve raw voltage as `V_a`;
- save diagnostic tables for point counts and window ranges.

### Step 4: Build static template library

Use the current proposed method's template-library functions from
`main_inv_log_2_method`:

- `load_inv_log_2_project_context.m`
- `make_response_template_library.m`
- `run_inv_log_2_main_method.m`

For the first feasibility check, use the existing unified static clearance
library. Do not attempt new static calibration yet.

### Step 5: Run three comparison methods

For each selected blade and time window, compare:

1. fixed-clearance template fit;
2. vibration-only fit with fixed clearance;
3. proposed joint clearance-vibration fit.

The minimum output per window should include:

- estimated clearance;
- estimated displacement offset;
- identified frequency or frequencies;
- identified amplitude;
- waveform RMSE;
- selected branch or candidate diagnostics;
- number of valid waveform samples.

### Step 6: Add independent frequency reference

For the same time window, extract the strain-gauge spectrum:

- STFT ridge or FFT peak;
- dominant frequency band;
- order relation with rotor speed.

The first validation metric is frequency agreement:

```text
|f_waveform - f_strain| / f_strain
```

### Step 7: Feasibility decision

The method is considered feasible on the existing experiments if:

- the adapter produces stable `highMap` inputs with enough waveform points;
- the joint model reduces waveform RMSE relative to the fixed template;
- identified frequency agrees with strain-gauge reference in the resonance
  window;
- estimated clearance does not drift wildly between adjacent windows unless
  the experiment intentionally changed the probe gap;
- residual plots do not show systematic missing waveform shape.

## First Minimal Validation Plan

### Priority 1: 20251222, blade 1, sensors 1-3

Suggested initial case:

- dynamic case: `1000_2500_3500`
- target blade: `1`
- sensors: `1, 2, 3`
- start time: `50 s`
- laps: `20`

Expected outputs:

- `outputs/20251222/B1_S123_highmap.mat`
- `outputs/20251222/B1_S123_result_joint.mat`
- `outputs/20251222/B1_S123_summary.csv`
- diagnostic figure comparing measured and fitted waveform samples.

### Priority 2: 20251222 sliding windows

After the first case runs, scan several windows around the suspected
resonance region and build a trend table:

- clearance estimate vs time;
- frequency estimate vs time;
- waveform RMSE vs time;
- strain-gauge reference frequency vs time.

### Priority 3: 20250527 capacitance sensors

Repeat the same workflow using only capacitance probes:

- dynamic case: `20250526_2500-3500_t400`
- sensors: `1, 3, 6`
- target blades: start with `1`, then extend to all six blades if stable.

## Figures To Produce

1. Real dynamic waveform construction:
   OPR, RPM, selected pulse windows, and mapped local waveforms.
2. Waveform fit comparison:
   measured waveform vs fixed template vs joint model.
3. Frequency reference comparison:
   identified waveform frequency vs strain-gauge spectrum.
4. Clearance-vibration crosstalk diagnostic:
   fixed-template apparent displacement drift vs joint-model clearance and
   vibration estimates.
5. Repeatability or window trend:
   estimates across adjacent windows or repeated blade passages.

## Interpretation Rules

- Do not claim absolute vibration amplitude truth unless an independent
  displacement reference is available.
- It is acceptable to claim frequency validation using strain-gauge data.
- It is acceptable to claim reduced clearance-vibration crosstalk if the
  joint model lowers residuals and stabilizes vibration estimates under
  clearance changes.
- If one sensor behaves differently, report it as sensor-specific bias or
  installation/gap difference rather than averaging it away.
- Keep 20250527 eddy-current channels separate from capacitance channels
  unless their response model is explicitly rebuilt.

## Current Folder Layout

```text
experimental_dynamic_decoupling_validation/
  README.md
  20251222/
    legacy/   preprocessing only
      output/ copied old preprocessing products
    outputs/
    Step01_Extract_OPR_Blade_Timing_20251222.m
    Step02_Calc_FullTime_BTT_Displacement_20251222.m
    Step03_Locate_Vibration_By_Strain_OPR_20251222.m
    Step04_Detect_Resonance_By_BTT_STE_20251222.m
    README_20251222.md
  20250527/
    legacy/   preprocessing only
      output/ copied old preprocessing products
    outputs/
    Step01_Extract_OPR_Blade_Timing_20250527.m
    Step02_Calc_FullTime_BTT_Displacement_20250527.m
    Step03_Locate_Vibration_By_Strain_OPR_20250527.m
    Step04_Detect_Resonance_By_BTT_STE_20250527.m
    README_20250527.md
```

## Programming Rule

All new validation programs should be written as independent Step scripts,
not as a short entry script plus many hidden subfunctions. Each Step script
should be readable from top to bottom with clear `%%` sections:

```matlab
%% 1. Paths and parameters
%% 2. Load existing data
%% 3. Main calculation
%% 4. Visualization
%% 5. Save outputs
```

Functions are allowed only for the copied, already-verified legacy
preprocessing programs or for truly shared utilities after the Step logic is
stable. This keeps the experimental workflow easy to inspect and debug.

## Current First-Batch Programs

The active first batch now follows the verified legacy preprocessing
programs and copied legacy preprocessing products for each experiment. The
copied scripts and `legacy/output/` files are treated as the trusted
data-reading, OPR, blade-indexing, waveform-windowing and strain-OPR
diagnostic layer.

Old Step4 static calibration and Step5 super-Gaussian identification are
intentionally excluded. Those stages will be rewritten for the proposed
clearance-vibration decoupling method.

Run from MATLAB:

```matlab
cd('E:\0小论文+程序\0博士期间小论文+程序\间隙和振动解耦\experimental_dynamic_decoupling_validation\20251222')
Step01_Extract_OPR_Blade_Timing_20251222
Step02_Calc_FullTime_BTT_Displacement_20251222
Step03_Locate_Vibration_By_Strain_OPR_20251222

cd('E:\0小论文+程序\0博士期间小论文+程序\间隙和振动解耦\experimental_dynamic_decoupling_validation\20250527')
Step01_Extract_OPR_Blade_Timing_20250527
Step02_Calc_FullTime_BTT_Displacement_20250527
Step03_Locate_Vibration_By_Strain_OPR_20250527
```

Expected products:

- `outputs/Step01_OPR_Blade_Timing_*.png`;
- `outputs/Step02_FullTime_BTT_Displacement_*.mat`;
- `outputs/Step02_FullTime_BTT_Displacement_*.png`;
- `outputs/Step03_Strain_OPR_Vibration_Window_*.mat`;
- `outputs/Step03_Strain_OPR_Vibration_Window_*.png`;
- `outputs/Step04_BTT_STE_Resonance_Regions_*.mat`;
- `outputs/Step04_BTT_STE_Resonance_Regions_*.csv`;
- `outputs/Step04_BTT_STE_Resonance_Regions_*.png`;
- `outputs/Step04_Aligned_Strain_STFT_Validation_*.png`.

These products define the candidate time windows that will be passed to the
next-stage BTT waveform extraction and decoupling analysis.

## Resonance Region Rule

Step04 determines candidate resonance regions by combining two signals:

1. BTT displacement short-time energy, after removing slow blade/sensor
   bias with a moving median;
2. strain-gauge synchronous-order energy from Step03, projected onto the
   BTT time axis.

A region is accepted only when the combined score is high and the
synchronous order frequency `best_order * RPM / 60` falls inside the
strain-observed frequency band, currently `450-700 Hz`. This avoids treating
large low-speed or run-down noise as resonance.
