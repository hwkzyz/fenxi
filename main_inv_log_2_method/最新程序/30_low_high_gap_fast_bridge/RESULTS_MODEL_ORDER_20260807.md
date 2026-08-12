# Adaptive single/dual model-order route (2026-08-07)

The formal entry now distinguishes model order before launching the expensive
dual-frequency search:

```text
independent single-frequency VP
-> full-voltage single refinement
-> projected residual second-frequency test
-> strong single: stop at single model
-> possible dual: run dual VP and full dual refinement
```

The residual test projects candidate second-frequency columns off gap,
translation, amplitude, phase, and first-frequency tangent directions. It uses
the complete trusted waveform and compares the conditional SSE reduction with a
three-parameter BIC penalty. A weak second component is not discarded by point
selection; it is reported through its conditional amplitude z-score.

## Regression

`Run_FormalMainRegression(true)` remains 6/6 successful. The three single cases
are classified as strong single and skip dual search. The three dual cases are
sent to the full dual route and remain successful.

At 5 dB, the off-grid dual regression remains 3/3 successful for both tested
seeds. All cases use 8912 trusted samples.

The amplitude audit in `Run_ModelOrderAmplitudeRegression` tested dual
amplitudes `0.15`, `0.08`, and `0.04 mm` against a `0.25 mm` primary component at
20 and 5 dB. All six cases entered the dual route and were recovered. The
smallest tested component (`0.04 mm`, 5 dB) had a conditional amplitude
z-score of 6.56 and a maximum frequency error of 0.294 Hz.

Typical single-case time fell to about 3.5--4.0 s in the full regression,
while dual cases remained about 14--19 s. The final full-voltage objective and
all trusted samples are unchanged.

The second-frequency test is deliberately conservative: only a positive
single-model evidence margin and a low conditional amplitude z-score allow the
dual search to be skipped. Ambiguous cases retain the dual safety search.
