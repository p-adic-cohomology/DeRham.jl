@testset "reduction engine interfaces" begin
    P = DeRham.HypersurfaceQuotient(2, 3)
    geometry = DeRham.ReductionEngineGeometry(P)
    @test geometry.reduction_degree == DeRham._reduction_g_degree(P)
    @test geometry.terminal_order == DeRham.terminal_pole_order(P)

    key = DeRham.ReductionChainKey([1, 2, 3], 4)
    @test key.u == (1, 2, 3)
    @test key.order == 4
    @test key == DeRham.ReductionChainKey([1, 2, 3], 4)

    handle = DeRham.ReductionChainHandle(key, 7)
    @test handle isa DeRham.AbstractReductionChainHandle
    @test handle.key === key
    @test handle.identifier == 7

    directions = DeRham.ReductionDirectionSet((u, k, v) -> sum(v) == 3 && all(u .>= k .* v))
    @test directions([3, 3, 3], 1, [1, 1, 1])
    @test !directions([3, 3, 3], 1, [2, 1, 1])

    @test isabstracttype(DeRham.AbstractReductionEngine)
    @test isabstracttype(DeRham.AbstractReductionPolicy)
    @test isdefined(DeRham, :materialize)
    @test parentmodule(DeRham.materialize) === DeRham
    @test isdefined(DeRham, :reduction_chain_key)
    @test isdefined(DeRham, :reduce!)
    @test isdefined(DeRham, :reduce_batched!)
    @test isdefined(DeRham, :combine!)
    @test isdefined(DeRham, :finish!)
    @test isdefined(DeRham, :run_reduction_policy!)
end
