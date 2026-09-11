# step07jcore MATLAB helpers for 20241106

This folder-local package makes the 20241106 gap-aware route self-contained.
It is copied here intentionally instead of being loaded from a root-level
shared package.

Current shared contract:

- Step07J mainline loads fixed joint static eta from the folder-local
  `Step07J_FixedJointEtaPrior_20241106.m`; MAT-file loading is retained
  only for explicit `STEP07J_STATIC_ETA_FILE` overrides;
- VP seeds use displacement-domain observations,
  `q_obs = -(V - F0)./Fx - eta_s`;
- weak template gradients are gated and down-weighted before VP fitting;
- when a full waveform evaluator is supplied, EO seeds are ranked by the
  complete weighted waveform RMSE, not by the linear VP residual.

The active callers are the three `Step07J_NestedStaticWarp_VPFullWave_*`
gap-aware scripts. The archived pre-sync copies are kept under
`_archive_code_snapshots/20260708_pre_foundation_gap_sync`.

After editing helper code, run
`experimental_dynamic_decoupling_validation/Audit_Step07JCore_Sync_AcrossGapFolders.m`.
