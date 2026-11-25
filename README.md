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
| `src/Table.jl` | bicubic Hermite tables of a model: `build_table`, `save_table`, `load_table`, `state`, `state_d`, `hydrogen` |
| `src/ExtendedEOS.jl` | a table with the model it was built from, for (T, P) outside it (`table_model`, `ExtendedEOS`, `extended_fields`) |
| `data/pf_Kurucz.input` | Kurucz partition functions |
| `python/witt.py`, `python/witt_bifrost.py` | Wittmann's Python EOS as distributed, and `witt_bifrost.witt()`: the same with the shim's patches (tolerances, the gasc H₂ fix), so Python and `libwitt` agree to the solver tolerance (from HSEMHS's oracle) |

## Build

`Pkg.build("WittBifrostEOS")` (or `WittBifrostEOS.build_libwitt()`) generates the shim and compiles `libwitt` with `CXX` (default
`g++`, `clang++` on macOS) into a scratch space keyed by the shim's checksum (`libwitt_build_dir()`), never into the package
directory, which Pkg installs read-only. Build on the node type that loads it. `build_libwitt(dir)` builds elsewhere; point
`WITTBIFROSTEOS_LIBWITT` or `use_libwitt!(path)` at the result. Every user in a process loads one library (`libwitt_path`,
`libwitt_handle`).

On hyades, build with plain `g++` (`CXX=g++`): the cluster's `loadgcc.sh` sets `CXX` to an MPI wrapper that links the
hcoll OpenMPI, and a `libwitt` built that way fails to load in processes that use another MPI.

## Using a table

```julia
using WittBifrostEOS
tab = load_table("eos_table_fullz.h5")             # HSEMHS's production table, or build_table(model, envelope)
X, Y = log(6000.0), log(1.0e4)                       # ln T [K], ln P [Pa]
lnrho(tab, X, Y)                                     # (ln ρ, status)
lnrho_dY(tab, X, Y)                                  # (ln ρ, ∂ln ρ/∂ln P, status)
state(tab, X, Y)                                     # (ln ρ, ln nₑ, ln χ₅₀₀, status)
state_d(tab, X, Y)                                   # the same values and their ∂/∂ln T, ∂/∂ln P (10 values with status)
```

Values are SI (ρ in kg m⁻³, nₑ in m⁻³, χ₅₀₀ in m⁻¹). Outside the table every call returns NaN with status
`EOS_OUT_OF_ENVELOPE`; `in_table` tests a point first. `state_d` returns the derivatives of the same bicubic Hermite
interpolant that `state` evaluates (continuous across cells), which is what an adjoint of a solver using the table needs:
HSEMHS's equilibrium adjoint uses it. At the table's slope breaks (node pairs 10⁻⁹ apart in ln T) the derivative is that of
the cell the point falls in. NLTE tables (`eos_model(tab) == "nlte"`) also give `hydrogen(tab, X, Y)`, (ln n_HI, ln n_p, status).

## Tests

`Pkg.test()`: the shim's patches, the model (LTE and with the NLTE correction), and an NLTE table built, saved and loaded,
with the extended EOS outside it and `state_d` against central differences of `state` (relative error below 10⁻⁶). The Python
port has no test here; HSEMHS's oracle tests use it.

## Users

HSEMHS and FFNOInversion pin this package by Git revision (`[sources]` in their `Project.toml`); a change here reaches them
when their pins are moved. As of 25 Nov 2026 both pin 8e5d56a (`state_d`), merged as #1.

## Origin

The Julia code and the shim come from HSEMHS (`HSEMHS/src/eos`, `HSEMHS/deps/witt`); the HSEMHS implementation document's
sections are cited in the comments. Moving them here changes no result: the model and both production tables evaluate bit for
bit as in HSEMHS on a (T, P) grid from 2300 K to 30 MK and 10⁻⁴ to 10⁶ Pa.
