# Stage timings for the small input; not a substitute for public-API timings.
using DeRham, Oscar
include("fixtures.jl")
f = dense_fixture(only(filter(c -> c.p == 11, DENSE_DWORK_FIXTURES)))
function stages(f)
    p = 11
    validation = @elapsed begin
        DeRham._dense_dwork_input(f)
        @assert DeRham.check_smoothness(f)
    end
    native, reconstruction, assembly = 0.0, 0.0, 0.0
    for N = 3:5
        native += @elapsed A = DeRham._dense_dwork_matrix(f, p, N)
        reconstruction +=
            @elapsed c = DeRham._dense_dwork_lift(DeRham._dense_dwork_residues(A, p, N), p)
        if c !== nothing
            assembly = @elapsed DeRham.zeta_function(ZZ.(reverse(c)), p, 2)
            return (; validation, native, reconstruction, assembly)
        end
    end
    error("ambiguous lift")
end
stages(f)
for _ = 1:3
    println(stages(f))
end
