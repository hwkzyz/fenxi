# Static gap waveform library for dynamic decoupling validation

This folder stages the multi-gap reference-blade waveform library used by
`reference_blade_gap_analysis` so that Step05 can build a response surface
inside `experimental_dynamic_decoupling_validation`.

## Original source

- Source CSV: `E:\试验数据\20260521标定数据\20260524标定_低采样率\Data\low_sample_gap_visualization\same_blade_waveforms_across_gaps_2000Hz.csv`
- Source analysis loader: `E:\0小论文+程序\0博士期间小论文+程序\间隙和振动解耦\reference_blade_gap_analysis\load_common_data.m`

## Gap definition

The source CSV stores the recorded mechanical gap as `gapMm`.  The paper
analysis applies a zero correction:

```text
trueReferenceGapMm = recordedGapMm - 0.2 mm
```

Therefore Step05 should use `trueReferenceGapMm` as the response-surface
gap coordinate.

## Contents

- `raw_source/same_blade_waveforms_across_gaps_2000Hz.csv`: untouched copy of
  the original multi-gap waveform table. This table was generated from
  `selectedChannel = 1` in the source MATLAB script; that is MATLAB's 1-based
  channel index and is labelled as CH0 in the original plots.

- `raw_binary_source/`: copied original ARTSCOPE binary files, header files and
  MATLAB scripts used to regenerate/audit the CSV library.  Step05 does not
  need to read these files during normal response-surface fitting.
- `by_recorded_gap/*/waveforms.csv`: one CSV per recorded gap, with an added
  `trueReferenceGapMm` column.
- `gap_waveform_manifest.csv`: compact manifest of gap values and file paths.
- `SOURCE_METADATA.json`: machine-readable provenance and gap definitions.

## Columns

The waveform CSV columns include:

- `sourceName`
- `gapMm`
- `trueReferenceGapMm`
- `bladeId`
- `eventNo`
- `timeRelative_s`
- `voltageMv`

All waveform voltages in the staged CSV are from channel 1 in MATLAB indexing
(`selectedChannel = 1`). File names such as `ACTS1000_data_2mm` indicate the
recorded gap is 2 mm; they do not indicate channel 2.

The time coordinate is the peak-centered local coordinate used by the
reference-blade gap analysis.  Step05 should rebuild the trusted spatial/time
window and baseline correction before fitting \(F(g,x)\).


## Raw binary archive

The `raw_binary_source/` folder contains a copied, read-only archive of the
original files in:

```text
E:\试验数据\20260521标定数据\20260524标定_低采样率\Data
```

It includes the `ACTS1000_data*.bin` files, their `*_header.txt` files, and
the MATLAB scripts that were used to parse and visualize the low-sampling-rate
gap waveforms.  These files are kept for provenance and regeneration only.  The
normal Step05 path should read `gap_waveform_manifest.csv` and the staged CSV
waveforms under `by_recorded_gap/`.
