# The optional backend is built once, never from a solve. No Python or OpenMP.
source = joinpath(@__DIR__, "native", "dense_frobenius.cpp")
destination = joinpath(@__DIR__, "usr", "bin")
mkpath(destination)
compiler = get(ENV, "CXX", "c++")
for n in (4, 5, 6)
    name = n == 4 ? "dense_dwork" : "dense_dwork_$n"
    target = joinpath(destination, name * (Sys.iswindows() ? ".exe" : ""))
    try
        mktempdir(destination) do temp
            binary = joinpath(temp, basename(target))
            run(`$compiler -O3 -std=c++17 -DNV=$n $source -o $binary`)
            mv(binary, target; force = true)
        end
    catch exception
        isfile(target) && rm(target)
        @warn "Optional dense Dwork backend for $n variables was not built. Install a C++17 compiler with 128-bit integer support, then run Pkg.build(\"DeRham\")." exception
    end
end
