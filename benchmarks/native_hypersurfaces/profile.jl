# Separate public-call costs on the small archived cubic threefold.
using DeRham, Oscar
include("fixtures.jl")
function profile()
    case=only(filter(c->c.id=="native_dwork_n5_p13", NATIVE_HYPERSURFACE_FIXTURES))
    f=native_fixture(case)
    DeRham.dense_dwork_zeta(f)
    times=Dict(
        name=>Float64[] for
        name in (:validation, :smoothness, :matrix, :residues, :lift, :assembly)
    )
    for _ = 1:7
        GC.gc()
        push!(times[:validation], @elapsed p=DeRham._dense_dwork_input(f))
        push!(times[:smoothness], @elapsed @assert DeRham.check_smoothness(f))
        push!(times[:matrix], @elapsed A=DeRham._dense_dwork_matrix(f, p, 5))
        push!(times[:residues], @elapsed residues=DeRham._dense_dwork_residues(A, p, 5))
        push!(times[:lift], @elapsed c=DeRham._dense_dwork_lift(residues, p; dimension = 3))
        push!(times[:assembly], @elapsed Z=DeRham.zeta_function(ZZ.(reverse(c)), p, 3))
        @assert c==case.expected
    end
    for name in (:validation, :smoothness, :matrix, :residues, :lift, :assembly)
        println(name, " median_seconds=", sort(times[name])[4], " samples=", times[name])
    end
end
profile()
