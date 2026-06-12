# AGENTS.md — Project instructions for AI agents

## Environment activation

**Always activate the environment before running any commands:**

```sh
micromamba activate fortran
```

All build tools (`fpm`, `meson`, `gfortran`, etc.) are available inside this environment.

## Project overview

`hompack` is a Fortran library for finding zeros or fixed points of nonlinear systems using globally convergent probability-one homotopy algorithms. It is a modernization of the legacy HOMPACK90 Fortran 77/90 code.

- Language: **Fortran** (fixed-form `.f` sources being migrated to free-form `.f90`)
- Primary build system: **fpm** (Fortran Package Manager)
- Secondary build system: **meson**
- External dependencies: BLAS / LAPACK (linked via OpenBLAS by default)

## Source layout

```
src/          — library source files (.f and .f90)
test/         — test programs (TESTF.f, TESTP.f, TESTS.f)
example/      — example programs
external/     — bundled reference BLAS/LAPACK subset (fallback)
c/            — C API
```

## Building and testing

### With fpm (preferred)

```sh
micromamba activate fortran
fpm build --profile debug
fpm test --profile debug
fpm build --profile release
fpm test --profile release
```

Run a specific example:

```sh
fpm run --example "example_name" --profile debug
```

### With meson

```sh
micromamba activate fortran
meson setup builddir -Dbuild_tests=true
meson compile -C builddir
meson test -C builddir
```

Optional meson flags:

| Option           | Values                      | Default  |
|------------------|-----------------------------|----------|
| `blas_option`    | `auto`, `openblas`, `refblas` | `auto` |
| `build_tests`    | `true`, `false`             | `false`  |
| `build_examples` | `true`, `false`             | `false`  |
| `build_shared`   | `true`, `false`             | `false`  |

## Code conventions

- **Fortran standard**: Fortran 2018 target for new/modified code.
- **Source form**: Free-form (`.f90`) for new files; legacy fixed-form (`.f`) files are being converted.
- **Case**: Lower case for new/modified Fortran code.
- **Implicit typing**: Allowed in legacy files (`HOMPACK90.f`); avoid in new `.f90` files.
- **Interfaces**: Explicit interfaces required for BLAS and LAPACK calls — see `src/blas_interfaces.f90` and `src/lapack_interfaces.f90`.
- **Precision**: Use kinds defined in `src/hompack_kinds.f90` instead of hard-coded precision specifiers.
- The project links against **OpenBLAS** (`link = ["openblas"]` in `fpm.toml`).

## Important notes

- Do **not** edit files under `external/` — these are vendored reference BLAS/LAPACK routines.
- Do **not** edit files under `_site/`, `build/`, or `builddir/` — these are generated artifacts.
- The file `src/HOMPACK90.f` is the original legacy source kept for reference while migration is in progress.
- Tests use data files `test/INNHP.DAT` and `test/OUTHP.DAT`; do not delete them.

## Token optimization

Use these defaults unless the user explicitly asks for deep explanation:

- Keep responses short and outcome-first. Prefer concise bullet points over long narrative text.
- Minimize repeated summaries. Report only deltas after each change batch.
- Read files in targeted ranges (relevant functions/sections) instead of full-file reads when possible.
- Prefer fast, narrow searches (`rg` with specific patterns and include globs) before semantic exploration.
- Batch independent read-only tool calls in parallel to reduce turn count.
- Avoid re-reading unchanged files; cache context from earlier steps in the same task.
- Edit only the smallest required region; avoid unrelated refactors and formatting churn.
- Run focused verification first (single test/file/target), then broaden only if needed.
- When diagnostics are long, report top actionable errors first and omit repetitive lines.
- Do not include large command outputs verbatim; summarize key lines and conclusions.
- For plans, use short actionable steps and update status briefly.
- If uncertainty is low, act directly; ask clarifying questions only when blocked by ambiguity.
