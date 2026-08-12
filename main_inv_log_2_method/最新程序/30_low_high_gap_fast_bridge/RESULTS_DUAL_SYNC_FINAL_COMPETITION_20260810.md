# Final EO Competition Audit (2026-08-10)

## Purpose

Funnel V2 has already removed the candidate-pruning failure: every tested
pressure case retains the true EO pair through the first joint-GN stage and,
except for a small number at 5 dB, into the final three full refinements. This
audit asks a different question: when the true pair is evaluated by the final
solver, can it still lose to a competing pair on the noisy voltage waveform?

The three deliberately difficult cases were tested with 20 seeds at complete
DAQ-record SNR levels of 25, 15, and 5 dB:

- `C02_eo_10_12`: close pair `(10,12)`, amplitudes `(0.25,0.15) mm`.
- `X01_corner_near`: close pair `(10,12)`, amplitudes `(0.25,0.10) mm`,
  `0.5 -> 0.7 mm` gap transition.
- `X02_corner_lowamp`: pair `(15,18)`, amplitudes `(0.15,0.10) mm`,
  `0.7 -> 0.5 mm` gap transition.

Noise is added to the complete raw DAQ record before spatial mapping. Each
condition now has a stable numeric noise seed, so the same condition and seed
produce the same waveform whether it is run alone or inside a full matrix.

## Results

| DAQ SNR | C02 final EO | X01 final EO | X02 final EO | True pair in final Top-3 |
|---:|---:|---:|---:|---:|
| 25 dB | 20/20 | 0/20 | 20/20 | 100%, 100%, 100% |
| 15 dB | 15/20 | 1/20 | 9/20 | 100%, 100%, 100% |
| 5 dB | 11/20 | 9/20 | 3/20 | 90%, 100%, 95% |

The detailed numerical table is
`output/dual_sync_funnel_v2_generalization_daq_snr/summary_final_competition_v2.csv`.
The corresponding plot is `final_competition_v2.png`.

## Interpretation

1. The final errors are not caused by `K2=24`, Top-K pruning, or missing
   full-model evaluation. The true EO pair is already in the final set.
2. `X01_corner_near` is a genuine order-selection ambiguity for the current
   8-revolution, full-waveform measurement geometry. A competing pair has a
   slightly smaller noisy-data residual even at 25 dB.
3. The current residual-margin heuristic is not a valid acceptance rule. Its
   noise-normalized margin is below 2 for every one of the 180 pressure
   realizations, including all 40 correct 25 dB C02/X02 solutions. It must not
   be used to declare the selected EO pair "identified" or "ambiguous" in the
   formal method.

## Oracle-order control

For `X01_corner_near`, the solver was rerun with the true `(10,12)` pair as
the only permitted order model. Parameter recovery then succeeds in:

| DAQ SNR | Parameter success | Median maximum amplitude error |
|---:|---:|---:|
| 25 dB | 20/20 | 0.0026 mm |
| 15 dB | 17/20 | 0.0061 mm |
| 5 dB | 10/20 | 0.0110 mm |

Thus the forward model and continuous parameter optimization are sound for
this case when its EO structure is known. The failure is specifically the
selection between close, competing EO structures.

## Consequence for the frozen method

Keep Funnel V2 with `K2=24` as the formal efficient solver for the validated
regular operating range. Do not increase Top-K further to address these cases.
For unknown close-order records, the paper should require either an external
modal/order prior, a longer record, or additional measurement information; the
current 8-revolution single-waveform data cannot support a universal
close-order blind-identification claim.

The next useful study is therefore not another candidate-budget sweep. It is
an information study of record length, sensor count, and physically justified
EO priors for the close-order boundary.
