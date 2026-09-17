# Full-point efficiency update (2026-08-07)

## Constraint

The formal solver keeps every sample that passes the trusted spatial-domain
check. Information scores are diagnostic and do not reduce the waveform used
by VP screening, nonlinear refinement, or final voltage-RMSE ranking.

## Exact computational changes

1. Dual-frequency VP forms sufficient statistics from the complete waveform
   and solves frequency-pair normal equations in batches.
2. The independent single-frequency VP caches its Fourier basis and solves all
   3-by-3 frequency systems in pages instead of rebuilding a weighted design
   matrix for every gap-frequency pair.
3. Conditional dual completion caches the full-band Fourier basis and the
   forward template for each tested gap.
4. Static gap/dx grids evaluate all dx values as a matrix. The coarse stage now
   uses all trusted samples by default instead of a 4000-point subset.

These changes alter computation order only. They do not alter the objective,
candidate grid, trusted samples, or final ranking metric.

## Regression results

`Run_FormalMainRegression(true)` passed 6/6 cases using 8912 high-speed samples
per case:

| Case | Maximum frequency error (Hz) | Elapsed (s) |
|---|---:|---:|
| 531 | 0.0252 | 12.34 |
| 733 | 0.0001 | 12.82 |
| 917 | 0.0110 | 14.79 |
| 500+1300 | 0.0061 | 14.84 |
| 700+1200 | 0.0549 | 14.67 |
| 900+1400 | 0.0931 | 19.71 |

The 5 dB off-grid seed-2 regression also passed 3/3 cases with all 8912
samples. Its maximum frequency error was 0.151 Hz and elapsed time was
12.77--18.53 s per case.

In the matched 5 dB smoke comparison, the original all-point scalar solver
required 19.84 s (single) and 24.96 s (dual). The current exact batched solver
required 8.39 s and 12.74 s in the same two cases. Frequency estimates and
voltage RMSE were unchanged to the reported precision.

## Remaining boundary

The default absolute-template formal route is covered by these regressions.
The optional `low_increment` route still rejects gap candidates whose queried
physical gap becomes non-positive; that route requires a separate physical
gap-bound correction and is not part of this efficiency result.
