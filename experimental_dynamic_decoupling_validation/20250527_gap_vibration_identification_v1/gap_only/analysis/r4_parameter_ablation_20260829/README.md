# R4 parameter ablation on the 20250527 experiment

This folder is an analysis fork. It does not modify the formal `results` files or rerun `Main09`/`Main10`.

## What the existing experiment can test immediately

The formal gap-aware result already contains the paired `fixed` and `gap_only` fits on the same high-speed windows. This directly tests the effect of removing the static gap increment `deltaGapMm` while keeping the frozen low-speed template and the vibration model unchanged.

The existing result does **not** expose a free voltage gain `a` or additive offset `b`. The voltage has already passed through the low-speed template and baseline preparation. Therefore `a` and `b` cannot be ablated honestly from the saved high-speed result alone.

## Parameter interpretation for this experiment

- `b`: remove from the high-speed incremental model. It is a baseline term and cancels in `V(g+dg)-V(g)`. Re-estimating it in the high-speed fit would be an incorrect extra degree of freedom.
- `a`: keep as a low-speed inner scale unless the raw low-speed and high-speed voltage chain has an independent gain calibration or a documented common normalization. Setting `a=1` is a separate calibration assumption, not a harmless deletion.
- `g0`, `mu`, `tau`, `zeta`: freeze after low-speed calibration; their ablation must be performed by rebuilding the low-speed template/model, not by editing a high-speed fit result.
- `deltaGapMm`: this is the directly testable high-speed static-gap term in the existing formal output.

Run `Run_R4_ExperimentalAblationInventory` first. It inventories the formal fit fields and compares `fixed` versus `gap_only` on identical windows. A second stage should rerun the low-speed preparation with `a=1` and with the baseline removed, then feed both frozen models into the same `Main10` high-speed windows.
