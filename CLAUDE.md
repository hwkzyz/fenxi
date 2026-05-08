# Project Coding Guidelines

Use the following rules by default when writing, reviewing, or refactoring code in this project.

These guidelines bias toward caution over speed. For trivial edits, use judgment.

## 1. Think Before Coding

- Do not assume missing details silently.
- State assumptions explicitly when they matter.
- If multiple interpretations exist, surface them instead of picking one invisibly.
- If something is unclear, stop and ask.
- If a simpler approach exists, propose it.

## 2. Simplicity First

- Write the minimum code needed to solve the requested problem.
- Do not add features that were not requested.
- Do not introduce abstractions for one-off code.
- Do not add speculative flexibility or configurability.
- Avoid error handling for impossible scenarios.
- If a solution feels overengineered, simplify it.

## 3. Surgical Changes

- Touch only the code required for the task.
- Do not refactor unrelated code while implementing a request.
- Do not "clean up" adjacent code, comments, or formatting unless needed.
- Match the existing style of the codebase.
- Remove only the unused code created by your own changes.
- If unrelated dead code is noticed, mention it instead of deleting it.

## 4. Goal-Driven Execution

- Define a concrete success condition before implementing.
- Prefer verifiable outcomes over vague instructions.
- For bug fixes, reproduce the bug first when practical, then verify the fix.
- For refactors, preserve behavior and verify before and after.
- For multi-step tasks, make a short plan and verify each step.

## Default Working Style

- Prefer clear explanations over hidden assumptions.
- Prefer smaller diffs over broad rewrites.
- Prefer passing checks, tests, or other concrete validation over "should work".

## Project-Specific Rules

- This repository is primarily a research codebase for clearance-vibration decoupling, blade vibration identification, and related simulation and analysis workflows.
- MATLAB is the default implementation language for core research pipelines. Use Python mainly for support tooling, automation, data handling, or tests unless a Python rewrite is explicitly requested.
- Preserve the existing step-based pipeline structure. Files such as `Step_*`, `Run_All*`, and `Config` or `Step_0_Config*` should keep their orchestration roles clear.
- Prefer extending the nearest existing pipeline or method variant instead of creating a brand-new framework.
- Keep method assumptions explicit, especially for physical meaning, units, normalization choices, sensor geometry, static gap definitions, and vibration parameter definitions.
- When a numerical choice matters, document why it was chosen and how success will be checked.

## Experiment And Results Discipline

- Make experiment changes reproducible through code, not manual post-processing.
- Prefer writing outputs into the existing results folders for that pipeline instead of scattering files across the repository.
- Do not overwrite important result artifacts or baseline outputs unless the task clearly requires regeneration.
- When modifying an algorithm, verify with a representative scenario, metric, figure, or saved intermediate result whenever practical.
- For comparisons between methods, preserve fair inputs and comparable settings unless the task is specifically to change the comparison protocol.

## Safe Edit Boundaries

- Avoid editing archived or historical material such as `old_archive`, generated caches, vendored dependencies, or tool bootstrap folders unless the request is directly about them.
- Avoid broad renaming across step scripts unless needed for correctness.
- Keep `Run_All` entry points usable and easy to follow.
- If a script appears to be a one-off exploration or scratch file, prefer creating a clearly named new file only when necessary rather than silently repurposing it.

## MATLAB And Figure Conventions

- Keep MATLAB scripts readable and procedural unless a stronger abstraction is clearly justified.
- Prefer explicit variable names over compressed clever code in analysis steps.
- When generating figures for papers or reports, keep style, labels, legends, and exported outputs consistent with surrounding scripts and manuscript terminology.
- Do not change figure appearance purely for aesthetics if it may affect comparison with existing results.

## Paper And Terminology Alignment

- Match the terminology used in the manuscript and theory notes when editing code comments, figure labels, captions, or method descriptions.
- Distinguish clearly between validated conclusions, exploratory observations, and hypotheses.
- Do not strengthen claims beyond what the current experiments or derivations support.
