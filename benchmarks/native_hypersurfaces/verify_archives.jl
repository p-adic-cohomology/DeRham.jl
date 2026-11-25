using DeRham, Oscar, Test
include("fixtures.jl")
for case in NATIVE_HYPERSURFACE_FIXTURES
    f=native_fixture(case)
    seconds=@elapsed cs=DeRham.dense_dwork_coefficients(f)
    @test BigInt.(reverse(cs))==case.expected
    z=DeRham.dense_dwork_zeta(f)
    R, T=polynomial_ring(ZZ, "T")
    F=fraction_field(R)
    P=R(case.expected)
    ambient=prod(1-case.p^i*T for i = 0:(case.n-2))
    expected=isodd(case.n-2) ? F(P)/F(ambient) : F(1)/(F(P)*F(ambient))
    @test z==expected
    println(case.id, " all_coefficients=", length(cs), " seconds=", seconds)
    flush(stdout)
end
