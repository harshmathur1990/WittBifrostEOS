"""
    WittBifrostEOS

Wittmann's equation of state for solar atmospheres, shared by HSEMHS and FFNOInversion: the C++ code with HSEMHS's shim
(tolerances, the H₂ fix, a status-returning solve) and FFNOInversion's synthesis entry points (`libwitt_path`, `build_libwitt`);
the direct solver and robust hot branch (`EOSModel`, `model_fields`); the NLTE hydrogen ionization learned from Bifrost
non-equilibrium (hion) runs (`load_nlte_correction`, `default_nlte_correction`); χ₅₀₀; and bicubic Hermite tables of the model
(`build_table`, `save_table`, `load_table`, `state`, `hydrogen`), with the model outside them (`ExtendedEOS`).
"""
module WittBifrostEOS

using HDF5
using SHA
using Libdl
using Printf
using Random
using Statistics

include("Status.jl")
include("WittmannDirect.jl")
include("SahaHot.jl")
include("NLTE.jl")
include("Opacity500.jl")
include("Model.jl")
include("Table.jl")
include("ExtendedEOS.jl")

"""The Kurucz partition functions (`pf_Kurucz.input`) distributed with Wittmann's EOS."""
default_pf_file() = joinpath(pkgdir(@__MODULE__), "data", "pf_Kurucz.input")
"""The NLTE hydrogen correction from the Bifrost hion runs (scripts/make_nlte_correction.py)."""
default_nlte_correction() = joinpath(pkgdir(@__MODULE__), "data", "nlte", "nlte_hion_correction.h5")

export EOS_OK, EOS_HOT_BRANCH, EOS_BLEND, EOS_OUT_OF_ENVELOPE, OPACITY_INCOMPLETE, EOS_DIRECT_FAILED, EOS_ROBUST_FALLBACK, EOS_NLTE
export InvalidInput, K_B_WITT, M_U_WITT, meta_value
export WittCtx, WittState, witt_state, witt_partials!, witt_mass_per_H, witt_abundances, witt_pf, witt_fix_gasc_h2, file_sha256
export libwitt_path, libwitt_handle, build_libwitt, default_pf_file, default_nlte_correction
export HotEOS, NLTECorrection, load_nlte_correction, EOSModel, model_fields, model_eval, Envelope, coverage_scan, select_T_sw
export AbstractEOS, EOSTable, build_table, save_table, load_table, validate_table, derivative_check, table_sha256
export lnrho, lnrho_dY, state, hydrogen, fields_at, in_table, eos_model, mass_per_H
export read_log_abund, table_model, ExtendedEOS, extended_model, extended_fields

end
