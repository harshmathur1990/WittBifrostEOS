# Pkg.build hook: generate the shim and compile libwitt into deps/witt (deps/witt/build.jl).
include(joinpath(@__DIR__, "witt", "build.jl"))
build()
