# The optional backend is built once, never from a solve. No Python or OpenMP.
source = joinpath(@__DIR__, "native", "dense_frobenius.cpp")
destination = joinpath(@__DIR__, "usr", "bin")
mkpath(destination)
compiler = get(ENV, "CXX", "c++")
target = joinpath(destination, "dense_dwork" * (Sys.iswindows() ? ".exe" : ""))
try
    mktempdir(destination) do temp
        binary = joinpath(temp, basename(target))
        run(`$compiler -O3 -std=c++17 $source -o $binary`)
        mv(binary, target; force = true)
    end
catch exception
    # DeRham's existing Julia algorithms do not require a C++ toolchain.
    isfile(target) && rm(target)
    @warn "Optional dense Dwork backend was not built. Install a C++17 compiler with 128-bit integer support, then run Pkg.build(\"DeRham\")." exception
end
