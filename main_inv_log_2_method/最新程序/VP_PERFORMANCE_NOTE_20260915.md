# Unified blind backend: VP funnel and runtime audit

## Why the previous smoke run was slow

The smoke test contains one operating case only.  Its expensive part was not
the number of cases, but the candidate evaluation pattern:

1. Every `delta_g` and EO pair was converted to a displacement seed.
2. Every pair was immediately evaluated with the complete voltage waveform.
3. The fixed-path evaluator then performed two variable-gap response queries
   for every sample.  Without a cache these queries used a MATLAB loop.
4. The retained candidates were evaluated again by bounded nonlinear fitting.

Thus the scan multiplied the number of clearance/EO pairs by the number of
waveform samples before the Top-K reduction.

## Experimental VP principle

The experimental Step05 programs use the same separation of roles:

\[
q_i = -\frac{V_i-F_0(x_i)}{F_x(x_i)},\qquad
q_i \approx d + a\sin(EO\,\Theta_i)+b\cos(EO\,\Theta_i).
\]

For each discrete EO (and, in the gap-aware route, each discrete clearance
candidate), the coefficients `[d,a,b]` are obtained by weighted linear least
squares.  This is the variable-projection seed.  Only a small adaptive set of
the best seeds is sent to the complete direct-template voltage fit, where
`A=sqrt(a^2+b^2)`, phase, `dx`, frequency offset, and clearance are refined.

## Current unified simulation backend

`run_unified_blind_dynamic_backend.m` now follows that contract:

- the complete `delta_g x EO` grid performs only the displacement-domain VP;
- a wider VP funnel (`max(topK,5*topK)`) is retained;
- voltage seed RMSE is evaluated only for that funnel;
- only the final `topK` candidates enter bounded full-wave refinement;
- `eval_path_increment_template.m` uses the fixed path `(gap,x)` cache, with
  the old pointwise evaluator retained only as a fallback.

The smoke run is still one case (`66` gap candidates, `15` EO candidates,
`990` VP pairs, `topK=3`).  The measured backend timing is approximately:

| stage | time |
|---|---:|
| VP scan | 0.21 s |
| voltage seed funnel | included in total backend timing |
| Top-K nonlinear refinement | 0.46 s |
| backend total | 0.68 s |

MATLAB process startup and low-speed template construction are outside the
backend timing and account for most of the command-line wall time.

The current smoke output still recovers EO and clearance, but its amplitude
estimate is biased (`0.151` mm versus `0.200` mm).  This is an identification
accuracy issue to investigate separately; it is not a runtime issue and must
not be hidden by increasing the candidate search.
