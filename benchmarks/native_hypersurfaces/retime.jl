using DeRham, Oscar
include("fixtures.jl")
include("../dense_dwork/fixtures.jl")
# Prime the common Julia wrapper and Oscar code paths before recording runtimes.
DeRham.dense_dwork_zeta(dense_fixture(first(DENSE_DWORK_FIXTURES)))
for case in NATIVE_HYPERSURFACE_FIXTURES
    f=native_fixture(case)
    expected=DeRham.zeta_function(reverse(case.expected), case.p, case.n-2)
    elapsed=@elapsed @assert DeRham.dense_dwork_zeta(f)==expected
    println(case.id, " runtime_seconds=", elapsed)
    flush(stdout)
end
for case in DENSE_DWORK_FIXTURES
    f=dense_fixture(case)
    expected=DeRham.zeta_function(reverse(case.expected), case.p, 2)
    elapsed=@elapsed @assert DeRham.dense_dwork_zeta(f)==expected
    println("dense_dwork_p", case.p, " runtime_seconds=", elapsed)
    flush(stdout)
end
println(
    "julia_version=",
    VERSION,
    " oscar_version=",
    pkgversion(Oscar),
    " threads=",
    Threads.nthreads(),
)
