# julia --project=. benchmarks/dense_dwork/benchmark.jl native 101
using DeRham, Oscar
include("fixtures.jl")
mode = get(ARGS, 1, "native")
p = parse(Int, get(ARGS, 2, "101"))
case = only(filter(c -> c.p == p, DENSE_DWORK_FIXTURES))
f = dense_fixture(case)
solve =
    mode == "native" ? () -> DeRham.dense_dwork_zeta(f) :
    () -> DeRham.zeta_function(f; use_threads = false, use_gpu = false, verbose = 0)
expected = DeRham.zeta_function(reverse(case.expected), p, 2)
first = @elapsed @assert solve() == expected
println("first_solve_seconds=", first);
flush(stdout)
get(ENV, "DENSE_BENCH_COLD", "0") == "1" && exit()
times = Float64[]
for _ = 1:3
    GC.gc()
    push!(times, @elapsed @assert solve() == expected)
    println("warm_repetition_seconds=", times[end])
    flush(stdout)
end
println("warm_seconds=", times)
println("median_seconds=", sort(times)[2])
println("julia_version=", VERSION, " oscar_version=", pkgversion(Oscar))
println("threads=", Threads.nthreads(), " backend=", mode)
println("full_projective_zeta_equal=true")
