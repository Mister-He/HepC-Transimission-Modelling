# AGENTS.md

## Working rules

- Keep changes inside this repository unless the user explicitly requests an
  external worktree or other path.
- Unless the user clearly asks for them, create as few new files and folders as
  possible. Prefer updating an existing file or overwriting an existing result
  path.
- Every substantive repository directory must contain a concise Markdown guide
  describing its files, purpose, and generated-versus-source status. Before
  traversing a directory, read that guide first. Generated timestamped log
  directories are documented by their parent rather than receiving one guide
  per run.
- Do not delete or overwrite user changes without first identifying the exact
  scope. Keep generated results in existing result directories.

The model-specific contract and calibration rules are in
[`docs/AGENTS.md`](docs/AGENTS.md). The current task specification is in
[`docs/structured_prompt.md`](docs/structured_prompt.md).
