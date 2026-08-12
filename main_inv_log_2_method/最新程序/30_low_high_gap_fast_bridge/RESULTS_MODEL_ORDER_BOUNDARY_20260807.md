# Model-order and frequency boundary audit (2026-08-07)

All cases used the complete 8912-sample trusted waveform. No point reduction
was used.

## Boundary sweep

| Case | Estimated order | Maximum frequency error | Projected second-frequency evidence |
|---|---:|---:|---:|
| 700+720, 10 dB | 2 | 0.190 Hz | strong |
| 700+750, 10 dB | 2 | 0.159 Hz | strong |
| 700+800, 10 dB | 2 | 300 Hz | strong but wrong basin |
| 700+720, 5 dB | 2 | 0.130 Hz | strong |
| 700+750, 5 dB | 2 | 300 Hz | strong but wrong basin |
| 700+800, 5 dB | 2 | 0.138 Hz | strong |
| 700+1200, A2=0.02 mm, 5 dB | 1 | not recovered | ambiguous/weak |
| 700+1200, A2=0.02 mm, 0 dB | 2 | 500 Hz | noise valley |
| single 531, 5 dB | 1 | 0.069 Hz | no strong second component |
| single 531, 0 dB | 2 | false dual | noise overfit |

The failures are not monotonic in frequency separation. A 20 Hz pair can be
recovered while a 100 Hz pair can select a remote valley for one noise
realization. The relevant boundary is therefore the voltage-objective margin
between frequency basins, not separation alone.

## Oracle audit

| Case | Selected RMSE (V) | True-basin/oracle RMSE (V) | Interpretation |
|---|---:|---:|---|
| 700+800, 10 dB | 0.0091246 | 0.0091158 | true basin is better; candidate/refinement failure |
| 700+750, 5 dB | 0.016614 | 0.016620 | wrong basin is slightly better; objective ambiguity |
| weak 700+1200, 5 dB | 0.016072 | 0.016059 | true dual is slightly better; order/candidate threshold boundary |
| weak 700+1200, 0 dB | 0.028721 | 0.028726 | wrong basin is slightly better; not uniquely identifiable |
| single 531, 0 dB | 0.028396 | 0.028424 at truth | dual model overfits noise; order ambiguity |

The RMSE differences at the boundary are only several microvolts to tens of
microvolts. A minimum-RMSE selector cannot reliably recover the generating
frequency in these cases without additional observations or a defensible model
order/noise penalty.

## Current practical boundary

The tested route is reliable for the standard `0.25+0.15 mm` cases down to
5 dB and for a `0.04 mm` secondary component at 5 dB. At `0.02 mm` and 5 dB,
the second component is already on the order-selection boundary. At 0 dB,
both false-dual selection and remote dual-frequency valleys occur.

The next algorithmic correction should retain multiple independently ranked
single-frequency seeds before the second-frequency test, then add a calibrated
effective-sample-size or block-bootstrap order penalty. Candidate fixes can
recover cases where the oracle true basin is lower, but they cannot resolve
cases where a wrong basin has lower voltage RMSE than the truth.
