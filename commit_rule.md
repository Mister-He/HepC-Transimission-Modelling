# Commit Message Conventions

Use the following format:

> [!info] Commit form
> <type>(optional scope): <description>

| Prefix | Purpose | Example |
|---|---|---|
| `feat` | Add a new feature | `feat: add age-specific screening` |
| `fix` | Fix a bug | `fix: correct prevalence denominator` |
| `docs` | Update documentation | `docs: explain calibration assumptions` |
| `style` | Change formatting without changing code behavior | `style: align assignment operators` |
| `refactor` | Restructure code without changing behavior | `refactor: extract likelihood function` |
| `perf` | Improve runtime performance | `perf: vectorize population updates` |
| `test` | Add or update tests | `test: cover empty simulation output` |
| `build` | Update the build system or dependencies | `build: update Rcpp dependency` |
| `ci` | Update continuous integration configuration | `ci: run tests on GitHub Actions` |
| `chore` | Perform other maintenance tasks | `chore: update gitignore` |
| `revert` | Revert a previous change | `revert: remove age-specific screening introduced in abc1234` |


# Guidelines

- Follow existing repository conventions (if applicable) when they differ from these defaults.
- Keep each commit focused on one logically complete change.
- Write concise, specific descriptions using imperative verbs, such as "add", "fix", or "remove".
- Use an optional scope to identify the affected module, such as `fix(calibration): ...`.
- For complex changes, add a body after a blank line explaining why the change is needed.
- Mark breaking changes with `!` before the colon, such as `feat(api)!: ...`. Explain the incompatibility and required migration in a `BREAKING CHANGE:` footer.
- Do not label changes to model assumptions, scientific calculations, or intended outputs as `refactor`.
- Review the staged diff before committing and exclude unrelated changes, temporary files, and secrets.
