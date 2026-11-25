# mode: native, upstream-default, upstream-depthfirst; id is an archived fixture.
using DeRham, Oscar
include("fixtures.jl")
function benchmark(mode, id)
    case=only(filter(c->c.id==id, NATIVE_HYPERSURFACE_FIXTURES))
    f=native_fixture(case)
    expected=DeRham.zeta_function(reverse(case.expected), case.p, case.n-2)
    solve=mode=="native" ? ()->DeRham.dense_dwork_zeta(f) :
          ()->DeRham.zeta_function(
        f;
        S = collect(0:(case.d-1)),
        algorithm = mode=="upstream-depthfirst" ? :depthfirst : :default,
        use_threads = false,
        use_gpu = false,
        verbose = 0,
    )
    first=@elapsed @assert solve()==expected
    println("first_solve_seconds=", first)
    flush(stdout)
    get(ENV, "NATIVE_COLD", "0")=="1" && return
    times=Float64[]
    for i = 1:5
        GC.gc()
        push!(times, @elapsed @assert solve()==expected)
        println("warm_repetition_seconds=", times[end])
        flush(stdout)
    end
    println("median_seconds=", sort(times)[3])
    println(
        "julia_version=",
        VERSION,
        " oscar_version=",
        pkgversion(Oscar),
        " threads=",
        Threads.nthreads(),
    )
    println("full_zeta_equal=true mode=", mode, " case=", id)
end
benchmark(ARGS...)
