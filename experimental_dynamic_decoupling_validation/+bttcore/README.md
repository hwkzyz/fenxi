# bttcore shared MATLAB helpers

This package is now a reference copy only. The active gap-aware scripts use
their own folder-local `+step07jcore` package so each dataset folder remains
self-contained.

Current shared contract:

- fixed joint static eta is loaded from
  `output/step05_joint_static_eta_preview/JointStaticEtaPreview_*.mat`;
- VP seeds use displacement-domain observations,
  `q_obs = -(V - F0)./Fx - eta_s`;
- weak template gradients are gated and down-weighted before VP fitting;
- when a full waveform evaluator is supplied, EO seeds are ranked by the
  complete weighted waveform RMSE, not by the linear VP residual.

The archived pre-sync copies are kept under
`_archive_code_snapshots/20260708_pre_foundation_gap_sync`.
