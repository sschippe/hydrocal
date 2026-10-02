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
| `USE_FLINT` | `ON` | use FLINT for the exact 1F1/2F1 engines when found |
| `USE_STATIC_OPENMP` | `OFF` | link libgomp statically (needed for a DLL-free Windows build) |

### Windows (cross-compiling from Linux)

Install mingw-w64 (`apt install g++-mingw-w64-x86-64`), then:

```bash
cmake -B build-win -S . \
      -DCMAKE_TOOLCHAIN_FILE=cmake/toolchain-mingw-x86_64.cmake \
      -DCMAKE_BUILD_TYPE=Release
cmake --build build-win --target win-zip
```

This produces `build-win/hydrocal.exe` and assembles `hydrocal_win.zip`
in the source root. The executable is statically linked, so it imports
only `KERNEL32.dll` and `msvcrt.dll` and needs no runtime DLLs.

Targets:

| Target | Meaning |
|---|---|
| `win-project` | write `hydrocal.vcxproj`, `hydrocal.sln` and `ReadMe.txt` to `build-win/win-project/` |
| `win-zip` | package the executable and project files as `hydrocal_win.zip` |

`win-project` also runs on a native build, so the Visual Studio files can be
regenerated at any time. They are derived from the same `HYDROCAL_SOURCES`
list CMake builds, so they cannot drift out of sync. Note that
`hydrocal_win.zip` is deliberately *not* tracked by git: it is published as a
GitHub Release asset.

Boost and FLINT are normally unavailable for the mingw target, so the
matching features fall back to the portable implementations and
`WARNINGS_AS_ERRORS` is disabled for this build.

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

## Citation

If you use hydrocal in your research, please cite the software itself:

```bibtex
@software{Schippers2026hydrocal,
  author       = {Schippers, Stefan},
  title        = {{hydrocal}: hydrogenic atomic structure calculations
                  and radiative transition rates},
  version      = {1.0.1},
  year         = {2026},
  month        = oct,
  publisher    = {GitHub},
  license      = {MIT},
  repository   = {github.com/sschippe/hydrocal},
  url          = {https://github.com/sschippe/hydrocal},
  orcid        = {0000-0002-6166-7138},
  note         = {Revision 2132; release \url{https://github.com/sschippe/hydrocal/releases/tag/v1.0.2}}
}
```

`@software` requires [biblatex](https://ctan.org/pkg/biblatex). For plain BibTeX use the `@misc` fallback in [`CITATION.bib`](CITATION.bib). Machine-readable metadata for GitHub's *Cite this repository* widget is in [`CITATION.cff`](CITATION.cff); both files also carry the method papers below.

There is no dedicated hydrocal software paper. The articles under [References](#references) describe the underlying physics, not the code, and should be cited only when referring to that physics.

## References

- S. Schippers et al., J. Phys. B **28** (1995) 3271 — [doi:10.1088/0953-4075/28/15/017](https://doi.org/10.1088/0953-4075/28/15/017)
- S. Schippers et al., ApJ **555** (2001) 1027 — [doi:10.1086/321512](https://doi.org/10.1086/321512)
- S. Schippers et al., A&A **421** (2004) 1185 — [doi:10.1051/0004-6361:20040380](https://doi.org/10.1051/0004-6361:20040380)
- S. Schippers, JQSRT **219** (2019) 33 — [doi:10.1016/j.jqsrt.2018.08.003](https://doi.org/10.1016/j.jqsrt.2018.08.003)
- S. Schippers et al., ChemPhysChem **24** (2023) e202300061 — [doi:10.1002/cphc.202300061](https://doi.org/10.1002/cphc.202300061)

## Author

Stefan Schippers — <stefan.schippers@uni-giessen.de>
Institute of Experimental Physics I, Justus-Liebig-Universität Giessen
[ORCID 0000-0002-6166-7138](https://orcid.org/0000-0002-6166-7138)
