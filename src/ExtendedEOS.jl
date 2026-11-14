# A table with the model it was built from: the table inside its envelope, the model outside it (e.g. a corona above the
# table's top temperature). The model is rebuilt from the table's metadata, with the checksums it records.

"""Log10 abundances of `eos.abundances`: "wittmann_default" (nothing) or a file of 99 values."""
function read_log_abund(spec::AbstractString)
    spec == "wittmann_default" && return nothing
    v = parse.(Float64, split(read(spec, String)))
    length(v) == 99 || invalid("abundance file $spec must contain 99 log10 fractions")
    return v
end

"""
    table_model(tab) -> EOSModel

The EOS model `tab` was built from: its partition functions, abundances, tolerances, hot-branch switch and, for an NLTE table,
its NLTE correction, from the metadata, with the partition-function and correction checksums verified.
"""
function table_model(tab::EOSTable)
    meta = tab.meta
    file_sha256(meta["pf_file"]) == meta["pf_sha256"] ||
        invalid("$(meta["pf_file"]) differs from the partition-function file the EOS table was built with")
    base = EOSModel(String(meta["pf_file"]); prec = Float64(meta["prec"]), maxit = Int(meta["maxit"]),
                    log_abund = read_log_abund(String(meta["abundances"])))
    nlte = nothing
    if eos_model(tab) == "nlte"
        path = String(meta["nlte_correction"])
        file_sha256(path) == meta["nlte_correction_sha256"] || invalid("$path differs from the NLTE correction the EOS table was built with")
        nlte = load_nlte_correction(path)
    end
    return with_nlte(with_switch(base, Float64(meta["T_sw"]), String(meta["hot_branch"])), nlte)
end

"""A table and, built on first use, the model it was built from (`table_model`)."""
mutable struct ExtendedEOS
    table::EOSTable
    model::Union{Nothing,EOSModel}
    lock::ReentrantLock
end
ExtendedEOS(tab::EOSTable) = ExtendedEOS(tab, nothing, ReentrantLock())

"""The model of an `ExtendedEOS`, built once (thread-safe)."""
function extended_model(e::ExtendedEOS)
    m = e.model
    m === nothing || return m
    lock(e.lock) do
        e.model === nothing && (e.model = table_model(e.table))
    end
    return e.model
end

"""
    extended_fields(e, X, Y) -> (lnrho, lnne, lnchi, lnnhi, lnnp)

At X = ln T, Y = ln P: the table inside its envelope (n_HI, n_p only from an NLTE table, NaN from an LTE one), the model outside;
NaNs where the model fails.
"""
function extended_fields(e::ExtendedEOS, X::Float64, Y::Float64)
    tab = e.table
    if in_table(tab, X, Y)
        p, q, t, s, hX, hY = locate(tab, X, Y)
        f(k) = hermite(tab, k, p, q, t, s, hX, hY, false)[1]
        return nfields(tab) >= 5 ? (f(1), f(2), f(3), f(4), f(5)) : (f(1), f(2), f(3), NaN, NaN)
    end
    v, _, ok = model_fields(extended_model(e), X, Y)
    return ok ? v : (NaN, NaN, NaN, NaN, NaN)
end
