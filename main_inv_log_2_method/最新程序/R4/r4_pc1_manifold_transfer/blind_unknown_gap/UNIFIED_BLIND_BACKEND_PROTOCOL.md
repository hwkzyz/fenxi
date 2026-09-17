# Unified blind gap--vibration backend

The R2/R3/R4 simulations use one dynamic identification backend.  Only the
static response adapter changes between routes:

- R2: direct same-geometry response surface;
- R3: direct gap-library response surface;
- R4: PC1/manifold-transferred response surface.

The forward model is `low_increment`:

\[
F_s(\Delta g,x)=T_s^{low}(x)+\kappa_s
\{R_s(x,g_{b,s}+\Delta g)-R_s(x,g_{b,s})\}.
\]

The estimator receives only a low-speed state, a dynamic observation bundle,
and response functions `evalF` and `evalFx`.  Truth EO, gap, amplitude, and
phase are used only by waveform generation and final error evaluation.

## Dynamic sequence

1. Generate a truth-independent candidate gap grid.
2. For each gap and all EO candidates, calculate the gradient displacement
   seed.
3. Solve the VP linear problem for static offset and sine/cosine coefficients.
4. Rank candidates by the gap-aware seed voltage RMSE and retain Top-K.
5. Refine every retained candidate with the complete voltage waveform model.
6. Rank only by data-only full-wave SSE.  The VP score is never the final
   decision metric.

The three scores are intentionally separate:

```text
vpDispRmseMm       displacement-domain VP score
seedVoltageRmseMv  cheap gap-aware voltage score
fullWaveSseMv2     final data-only objective
```

## Frequency model

The main simulation uses blind integer EO candidates and locally refines a
continuous frequency offset:

\[
u(t)=a\sin[EO\,\Theta(t)+2\pi\delta f(t-t_c)]
    +b\cos[EO\,\Theta(t)+2\pi\delta f(t-t_c)].
\]

The reported amplitude and phase are derived from `(a,b)`.  `deltaF` is not a
truth prior.

## Identification status

The backend reports `identified` only when the best candidate has a finite
second-candidate margin above the configured threshold.  Otherwise it reports
`ambiguous_joint`.  R2/R3/R4 batch results must preserve this status instead of
forcing every low-residual fit to be identified.

## Verification

The independent validation entry point is:

`Main_10_UnifiedBlindBackendValidation.m`

It writes results to `output/unified_backend_full` when run with the default
full four-case configuration.  The original `Main_08`/`Main_09` files and
their output folders are not overwritten.
