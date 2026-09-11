# Common R5 dynamic worker

`R5_Run_SensorConditionedWindowIdentification.m` is the shared dynamic
identification worker for the 20241106 and 20250527 migration packages.
Those package entry points add this directory to the MATLAB path before
calling the worker.

The worker contract is:

- the preparation stage freezes the sensor-conditioned latent, registration,
  gain, low-speed template, and reference gap;
- each high-speed window keeps the configured number of EO candidates;
- EO is selected from the voltage objective after candidate fitting;
- each candidate may refine a continuous frequency around its integer EO
  seed;
- the vibration phase uses the accumulated rotor angle plus a time-domain
  detuning term;
- the dynamic state contains only `EO`, `A`, `phase`, `dx`, and independent
  `dg_s`; strain is not part of the objective;
- `referenceEO` is recorded for audit only. Set `lockReferenceEO=true` only
  in an explicitly labelled diagnostic run.

The 20251222 `Main06` implementation remains the authoritative reference for
the candidate-search and continuous-frequency contract. The two migration
packages use this common worker because their static observation operators
are sensor-conditioned and dataset-specific.

Large raw-data and MAT result files are not required in this source-only
directory. Their locations and required fields are recorded by each package's
preparation manifest.
