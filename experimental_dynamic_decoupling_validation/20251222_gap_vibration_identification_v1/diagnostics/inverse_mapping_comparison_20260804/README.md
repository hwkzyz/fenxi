# 20251222 inverse-mapping comparison

This directory is an isolated diagnostic implementation. It does not modify
`Main10_GapAware_FullWave_Identification_20251222.m`, the frozen VP config, or
formal result files.

The first runnable check is `Diagnostic_InverseMapping_UnitSmoke_20251222.m`.
It validates the generic two-branch bounded inversion and EO seed fit using a
synthetic bell-shaped response. It is not an experimental result.

The experimental comparison must supply an evaluator copied or extracted from
Main10 exactly. The evaluator must return voltage, derivative, and query-safe
domain information for each sensor. A simplified simulation response is not an
acceptable substitute.

The current diagnostic keeps exactly the first three inverse-mapping EO
candidates. Candidate-rho is still recorded as a diagnostic but no longer
expands the full-wave optimization set to five candidates.
