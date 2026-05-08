# Project CLAUDE Template

Use this file as the default `CLAUDE.md` for a new project, then customize the placeholder sections for the specific repository.

## Core Working Principles

These rules apply by default when writing, reviewing, or refactoring code in this project.

- Do not assume missing details silently.
- State important assumptions explicitly.
- If multiple interpretations exist, surface them instead of choosing one invisibly.
- Prefer the simplest solution that satisfies the request.
- Do not add speculative abstractions, flexibility, or features.
- Touch only the code required for the task.
- Match the existing code style unless there is a clear reason not to.
- Define a concrete success condition before implementing.
- Verify results with tests, scripts, metrics, figures, or other project-appropriate checks.

## Project Summary

- Project purpose: [describe the engineering or research goal]
- Primary language(s): [Python / MATLAB / C++ / etc.]
- Secondary tooling: [automation, plotting, data processing, docs]
- Main deliverables: [software / experiments / figures / paper / report / dataset]

## Preferred Working Style

- Prefer smaller diffs over broad rewrites.
- Prefer reproducible changes through code rather than manual editing.
- Prefer extending nearby code or the nearest existing workflow instead of creating a parallel framework.
- If a simpler approach exists, propose it before implementing a more complex one.

## Repository Structure Rules

- Keep existing entry points easy to follow.
- Preserve the meaning of top-level scripts, pipeline steps, configs, and results folders.
- Prefer editing the nearest existing module or pipeline instead of creating a new one.
- Avoid broad renaming unless required for correctness.

Project-specific structure notes:

- Main entry points: [e.g. `Run_All.m`, `main.py`, `train.py`]
- Config files: [list or describe]
- Results/output folders: [list or describe]
- Scratch or archive folders to avoid by default: [list or describe]

## Language-Specific Rules

### MATLAB

- Keep scripts readable and procedural unless stronger abstraction is clearly justified.
- Keep units, physical assumptions, normalization choices, and variable meanings explicit.
- Prefer stable figure generation and saved intermediate artifacts for reproducibility.

### Python

- Use Python mainly for tooling, automation, tests, and data processing unless core algorithms are meant to live there.
- Prefer clear scripts and modules over clever indirection.
- Add tests when fixing bugs or changing behavior when practical.

### Other Languages

- Add project-specific rules here.

## Experiment And Validation Rules

- Changes to algorithms should be validated with a representative case whenever practical.
- Comparisons between methods should keep inputs and settings fair unless the task is to change the protocol.
- Do not overwrite important baseline outputs unless regeneration is part of the task.
- Save outputs in the established results location for that workflow.
- Distinguish clearly between validated conclusions, exploratory observations, and hypotheses.

Validation checklist for this project:

- Main correctness check: [test / script / metric / visual inspection / benchmark]
- Secondary checks: [list]
- Expensive checks to avoid unless needed: [list]

## Figures, Reports, And Terminology

- Keep figure labels, legends, captions, and exported filenames consistent with project terminology.
- Match manuscript or report terminology when editing comments, labels, or descriptions.
- Do not strengthen claims beyond what experiments, derivations, or tests support.

## Safe Edit Boundaries

Do not edit these areas unless the task directly requires it:

- Archived folders: [list]
- Vendored dependencies or generated caches: [list]
- Historical outputs or frozen baselines: [list]

If unrelated dead code or cleanup opportunities are noticed, mention them instead of deleting them by default.

## Multi-Step Task Pattern

For non-trivial tasks, work in this pattern:

1. Clarify assumptions and define success.
2. Inspect the nearest relevant files or pipeline.
3. Make the smallest reasonable implementation change.
4. Verify with the strongest practical check.
5. Report what changed, what was verified, and any remaining risk.

## Optional Project Notes

- Coding style preferences: [describe]
- Naming conventions: [describe]
- Data handling constraints: [describe]
- Performance constraints: [describe]
- Platform constraints: [Windows / Linux / cluster / MATLAB version / Python version]
