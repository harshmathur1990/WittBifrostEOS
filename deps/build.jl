# Pkg.build hook: generate the shim and compile libwitt into the package's scratch space (WittBifrostEOS.libwitt_build_dir),
# never into the package directory, which Pkg installs read-only.
import WittBifrostEOS
WittBifrostEOS.build_libwitt()
