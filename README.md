# hydrocal

A collection of subroutines performing **hydrogenic atomic structure calculations**, such as
bound–bound and bound–free radiative transition rates, which are then used to calculate cross
sections and rate coefficients for radiative recombination. Convolution of cross sections into
rate coefficients is supported as well.

Development is driven by data-evaluation needs of storage-ring electron–ion collision
experiments.

## Requirements

| | |
|---|---|
| Compiler | C++17 |
| Build system | CMake ≥ 3.16 |
| Required | a C++17 compiler and `libm` |

All external dependencies are **optional** and enable extra functionality when found:

| Dependency | Enables |
|---|---|
| OpenMP | multi-threaded evaluation |
| Boost (headers) | Boost-based routines |
| FLINT | exact 1F1 and 2F1 summation engines |
| FLTK + Cairo | the `hydrocalGUI` graphical front end |

## Building

```bash
cmake -B build -S .
cmake --build build -j
```

The CLI executable is written to `build/hydrocal`.

CMake options:

| Option | Default | Meaning |
|---|---|---|
| `ENABLE_WARNINGS` | `ON` | common warning set (`-Wall -Wextra -Wshadow`, `/W4`) |
| `WARNINGS_AS_ERRORS` | `ON` | promote warnings to errors (`-Werror`, `/WX`) |

## Running

`hydrocal` reads its input on stdin and writes results to stdout:

```bash
./build/hydrocal < doc/examples/U92RR.hcin
```

Each test case in `doc/examples/` is a `.hcin` input file.

## Tests

16 CTest cases compare output against the reference data in `doc/examples/`:

```bash
ctest --test-dir build --output-on-failure
```

## Version numbering

Every binary embeds a revision number, printed in the start-up banner and in
the header of most output files:

```
 *  hydrocal, revision 2123                                   *
```

It is a monotonically increasing integer that **continues the Subversion
revision numbering**, so numbers stay unique and correctly ordered across the
SVN → git transition. The last Subversion revision of the imported tree was
**r2122**, and the counter carries on from there.

The counter is the single line in the **`REVISION`** file at the top of the
source tree. It is bumped automatically by the `pre-commit` hook in
`.githooks/`, which writes the new value and stages it as part of the same
commit — so the number a build reports always belongs to the tree it was
built from, and there is nothing to bump by hand.

Git does not track hook configuration, so activate it once per clone:

```bash
git config core.hooksPath .githooks    # or: cmake --build build --target setup-hooks
```

To set the number deliberately (e.g. to mark a release), edit `REVISION`
before committing; the hook increments from whatever is in the file.

Commits made with `--no-verify` skip the bump — increment `REVISION` by hand
in that case so the sequence has no gaps. Note also that the counter follows
the file, not the history, so it stays put across `rebase`, `reset` and
history rewriting — already-published numbers remain valid.

If the `REVISION` file is missing (a partial copy, for instance) the build
falls back to reconstructing the value as
`REVISION_BASE + (commits on the main line) - 1`, where `REVISION_BASE` is set
in `CMakeLists.txt`.

A second macro, `HYDROCAL_GITID`, records `git describe` output (tag, distance
and abbreviated hash, plus `-dirty` for uncommitted changes) for pinning an
exact commit in a bug report.

## Layout

| Path | Contents |
|---|---|
| `src/` | the hydrocal library and CLI/GUI sources |
| `doc/` | Doxygen sources (`doc/src`), example inputs and references (`doc/examples`) |
| `autostructure/` | the autostructure (AUTO) structure-generation codes and examples |
| `JAC/` | JAC relativistic atomic structure code and install helpers |
| `REVISION` | current revision number, bumped by the `.githooks/` commit hook |

## License

MIT — see [LICENSE](LICENSE).

## References

- S. Schippers et al., J. Phys. B **28** (1995) 3271 — [doi:10.1088/0953-4075/28/15/017](https://doi.org/10.1088/0953-4075/28/15/017)
- S. Schippers et al., ApJ **555** (2001) 1027 — [doi:10.1086/321512](https://doi.org/10.1086/321512)
- S. Schippers et al., A&A **421** (2004) 1185 — [doi:10.1051/0004-6361:20040380](https://doi.org/10.1051/0004-6361:20040380)
- S. Schippers, JQSRT **219** (2019) 33 — [doi:10.1016/j.jqsrt.2018.08.003](https://doi.org/10.1016/j.jqsrt.2018.08.003)

## Author

Stefan Schippers — <stefan.schippers@uni-giessen.de>
