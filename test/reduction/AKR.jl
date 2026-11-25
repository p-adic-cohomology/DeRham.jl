@testset "AKR algorithm selector support" begin
    _, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    hypersurface = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    for algorithm in (:akr, :akr_naive, :akr_levelbylevel)
        @test DeRham._supports_reduction_algorithm(hypersurface, algorithm)
    end
end

@testset "AKR plane curve full-precision equality" begin
    curve = DeRham.fermat_hypersurface(3, 3, 2)
    _, naive = DeRham.zeta_coefficients(curve; givefrobmat = true, algorithm = :akr_naive)
    @test naive == BigInt[2, 0, 1]
    @test DeRham.zeta_coefficients(curve; algorithm = :akr) == naive
    @test DeRham.zeta_coefficients(curve; algorithm = :akr_levelbylevel) == naive
    @test DeRham.zeta_coefficients(curve) == naive
    @test DeRham.zeta_coefficients(curve; algorithm = :akr_naive, termorder = :lex) == naive
    @test DeRham.zeta_coefficients(
        curve;
        algorithm = :akr_levelbylevel,
        termorder = :lex,
    ) == naive
end

@testset "AKR quartic K3 limited-precision matrices" begin
    _, (x, y, z, w) = polynomial_ring(GF(2), ["x", "y", "z", "w"])
    k3 = x^3*y + y^3*z + z^3*w + w^3*x + x*y*z*w
    default_matrix = DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = :default)
    naive_matrix = DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = :akr_naive)
    levelbylevel_matrix =
        DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = :akr_levelbylevel)
    @test naive_matrix == default_matrix
    @test levelbylevel_matrix == default_matrix
end

@testset "handoff order is bounded before cache construction" begin
    curve = DeRham.fermat_hypersurface(3, 3, 2)
    @test_throws ArgumentError DeRham.zeta_coefficients(
        curve;
        algorithm = :depthfirst,
        dense_reduction = :akr,
        handoff_pole_order = typemax(Int),
    )
end

@testset "AKR dense reduction split equality" begin
    _, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    f = x^3 + y^3 + z^3
    n = 2
    handoff_order = n + 1
    d = 3
    S = collect(0:n)
    params = DeRham.default_params()
    quotient = DeRham.HypersurfaceQuotient(n, d)
    cache = DeRham.controlled_reduction_cache(n, d, S, params, quotient)
    (_, f, pseudo_inverse) =
        DeRham.find_Ssmooth_model(f, quotient, 1, S, params, true, cache)

    precisionring, _ = residue_ring(ZZ, 5)
    precisionringpoly, precisionvars = polynomial_ring(precisionring, ["x", "y", "z"])
    f_lift = DeRham.liftCoefficients(precisionring, precisionringpoly, f)
    pseudo_inverse_space =
        matrix_space(precisionring, nrows(pseudo_inverse), ncols(pseudo_inverse))
    pseudo_inverse_lifted = pseudo_inverse_space()
    for i = 1:nrows(pseudo_inverse), j = 1:ncols(pseudo_inverse)
        pseudo_inverse_lifted[i, j] = ZZ(pseudo_inverse[i, j])
    end

    picache = DeRham._build_akr_picache(
        4,
        n,
        d,
        S,
        f_lift,
        pseudo_inverse_lifted,
        cache,
        params,
        5;
        target_order = n,
    )
    context = (; S, f = f_lift, p = 5, picache, params, terminal_order = n)
    x_lift = precisionvars[1]
    input = [[x_lift^9, 4]]

    @testset "one-shot and split reductions agree" begin
        for reduction in
            (DeRham.AKRNaiveDenseReduction(), DeRham.AKRLevelByLevelDenseReduction())
            direct = DeRham.dense_reduce(reduction, input, 4, 2, context)
            middle = DeRham.dense_reduce(reduction, input, 4, handoff_order, context)
            split = DeRham.dense_reduce(reduction, middle, handoff_order, n, context)
            @test direct == split
            @test only(direct)[2] == n
        end
    end
end

@testset "dense-reduction setup options" begin
    curve = DeRham.fermat_hypersurface(3, 3, 2)
    frontend = DeRham.ProjectiveHypersurface(curve)
    n = 2
    handoff_order = n + 1
    default = DeRham.reduction_setup(frontend, [-1], :default)
    preset = DeRham.reduction_setup(frontend, [-1], :akr)
    mixed = DeRham.reduction_setup(
        frontend,
        [-1],
        :depthfirst;
        dense_reduction = :akr,
        handoff_pole_order = handoff_order,
    )

    @test default.dense_reduction === nothing
    @test default.handoff_pole_order === nothing
    @test preset.policy == :akr
    @test preset.dense_reduction == :akr_naive
    @test preset.handoff_pole_order === nothing
    @test mixed.policy == :depthfirst
    @test mixed.dense_reduction == :akr_naive
    @test mixed.handoff_pole_order == handoff_order
    @test_throws ArgumentError DeRham.reduction_setup(
        frontend,
        [-1],
        :depthfirst;
        handoff_pole_order = 3,
    )
    @test_throws ArgumentError DeRham.reduction_setup(
        frontend,
        [-1],
        :depthfirst;
        dense_reduction = :akr_levelbylevel,
        handoff_pole_order = 1,
    )
end
