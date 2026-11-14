# EOS status bits, the input error, and the constants the Wittmann C++ code uses (from HSEMHS, implementation doc §3.4, §4.5).

# per-node EOS bit mask
const EOS_OK              = 0x00
const EOS_HOT_BRANCH      = 0x01
const EOS_BLEND           = 0x02
const EOS_OUT_OF_ENVELOPE = 0x04
const OPACITY_INCOMPLETE  = 0x08
const EOS_DIRECT_FAILED   = 0x10

"""Error raised at API boundaries for rejected input. Carries a precise message."""
struct InvalidInput <: Exception
    msg::String
end
Base.showerror(io::IO, e::InvalidInput) = print(io, "invalid_input: ", e.msg)

invalid(msg) = throw(InvalidInput(msg))

# Constants used inside the Wittmann C++ code (CODATA 2010, CGS there), in SI.
# Quantities that must agree with the C++ EOS state use these.
const K_B_WITT = 1.3806488e-23
const M_U_WITT = 1.660538921e-27

# HDF5 metadata writer (tables; HSEMHS outputs)
meta_value(v::Symbol) = String(v)
meta_value(v::Nothing) = "none"
meta_value(v::AbstractVector{<:AbstractString}) = String.(v)
meta_value(v::Bool) = Int8(v)
meta_value(v) = v

"""Write a (nested) metadata Dict into HDF5 group `g`: dicts as groups, matrices and long vectors as datasets, the rest as attributes."""
function _write_meta(g, meta::AbstractDict)
    for (k, v) in meta
        k = String(k)
        if v isa AbstractDict
            sub = haskey(g, k) ? g[k] : create_group(g, k)
            _write_meta(sub, v)
        elseif v isa AbstractMatrix
            g[k] = permutedims(Matrix(v), (2, 1))
        elseif v isa AbstractVector && !(eltype(v) <: AbstractString) && length(v) > 64
            g[k] = collect(v)
        else
            write_attribute(g, k, meta_value(v))
        end
    end
end
