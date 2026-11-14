# WittBifrostEOS

Wittmann's equation of state for solar atmospheres, with the hydrogen ionization of Bifrost non-equilibrium (hion) runs as an
option. One implementation, used by HSEMHS (hydrostatic and magnetohydrostatic reconstruction) and FFNOInversion (synthesis and
inversion).

## Contents

| | |
| --- | --- |
| `deps/witt/witt_eos_cpp.cpp` | Wittmann's EOS (C++), as distributed, with FFNOInversion's synthesis entry points appended |
| `deps/witt/build.jl` | generates `witt_shim.cpp` by exact, counted patches (tolerances 1e-10, the gasc H₂ fix, a status-returning solve, the context C API) and compiles `libwitt` |
| `src/WittmannDirect.jl` | the ccall wrapper (`WittCtx`, `witt_state`, `witt_partials!`, …) |
| `src/SahaHot.jl` | the robust hot branch: the same model, solved without the C++ clamps |
| `src/NLTE.jl` | the NLTE hydrogen correction D(T, P) learned from Bifrost hion runs (`data/nlte`, `scripts/make_nlte_correction.py`) |
| `src/Opacity500.jl` | χ₅₀₀ |
| `src/Model.jl` | `EOSModel`, `model_fields`: direct solver, hot branch and their blend at T_sw; NLTE where D ≠ 0 |
| `src/Table.jl` | bicubic Hermite tables of a model: `build_table`, `save_table`, `load_table`, `state`, `hydrogen` |
| `src/ExtendedEOS.jl` | a table with the model it was built from, for (T, P) outside it (`table_model`, `ExtendedEOS`, `extended_fields`) |
| `data/pf_Kurucz.input` | Kurucz partition functions |

## Build

`Pkg.build("WittBifrostEOS")` (or `WittBifrostEOS.build_libwitt()`) compiles `deps/witt/libwitt.so` with `CXX` (default `g++`,
`clang++` on macOS). Where the package directory is read-only, build elsewhere with `build_libwitt(dir)` and set
`WITTBIFROSTEOS_LIBWITT` to the library. Every user in a process loads the same library (`libwitt_path`, `libwitt_handle`).

## Origin

The Julia code and the shim come from HSEMHS (`HSEMHS/src/eos`, `HSEMHS/deps/witt`); the HSEMHS implementation document's
sections are cited in the comments. Moving them here changes no result: the model and both production tables evaluate bit for
bit as in HSEMHS on a (T, P) grid from 2300 K to 30 MK and 10⁻⁴ to 10⁶ Pa.
