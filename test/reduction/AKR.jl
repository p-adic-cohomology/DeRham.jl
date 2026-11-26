@testset "AKR algorithm selector support" begin
    _, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    hypersurface = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    for algorithm in (:akr, :akr_naive, :akr_levelbylevel)
        @test DeRham._supports_reduction_algorithm(hypersurface, algorithm)
    end
end

@testset "AKR plane curve known output" begin
    curve = DeRham.fermat_hypersurface(3, 3, 7)
    expected_coefficients = BigInt[7, 1, 1]
    @test DeRham.zeta_coefficients(curve; algorithm = :akr) == expected_coefficients

    curve_p2 = DeRham.fermat_hypersurface(3, 3, 2)
    expected_lex_matrix = BigInt[0 1; 2 0]
    naive_lex_matrix = DeRham.frobenius_matrix_with_precision(
        curve_p2,
        1;
        algorithm = :akr_naive,
        termorder = :lex,
    )
    levelbylevel_lex_matrix = DeRham.frobenius_matrix_with_precision(
        curve_p2,
        1;
        algorithm = :akr_levelbylevel,
        termorder = :lex,
    )
    @test [
        BigInt(lift(naive_lex_matrix[i, j])) for
        i in axes(naive_lex_matrix, 1), j in axes(naive_lex_matrix, 2)
    ] == expected_lex_matrix
    @test [
        BigInt(lift(levelbylevel_lex_matrix[i, j])) for
        i in axes(levelbylevel_lex_matrix, 1), j in axes(levelbylevel_lex_matrix, 2)
    ] == expected_lex_matrix
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

@testset "handoff cache extends existing basis data" begin
    curve = DeRham.fermat_hypersurface(3, 3, 2)
    frontend = DeRham.ProjectiveHypersurface(curve)
    quotient = DeRham.cohomology_quotient(frontend)
    params = DeRham.default_params()
    handoff_order = DeRham.terminal_pole_order(quotient) + 1
    handoff_params = DeRham.ZetaFunctionParams(
        params.verbose,
        params.givefrobmat,
        :depthfirst,
        params.termorder,
        params.vars_reversed,
        params.fastevaluation,
        params.always_use_bigints,
        params.use_gpu,
        params.use_threads,
        :akr_naive,
        handoff_order,
    )
    S = collect(0:quotient.n)
    cache = DeRham.controlled_reduction_cache(quotient.n, quotient.d, S, params, quotient)
    terminal_degree =
        DeRham.numerator_degree(quotient, DeRham.terminal_pole_order(quotient))
    geometry = DeRham.ReductionEngineGeometry(quotient, handoff_params)
    handoff_degree = DeRham._reduction_g_degree(quotient, handoff_params)
    cached_terminal_basis = cache[terminal_degree]

    @test geometry.reduction_degree == DeRham._reduction_g_degree(quotient)
    @test geometry.reduction_degree != handoff_degree
    @test geometry.terminal_order == handoff_order

    @test DeRham.extend_controlled_reduction_cache!(
        cache,
        quotient.n,
        quotient.d,
        S,
        handoff_params,
        quotient,
    ) === cache
    @test cache[terminal_degree] === cached_terminal_basis
    @test handoff_degree in DeRham.cached_degrees(cache)

    expected_forward, expected_reverse = DeRham._controlled_reduction_cache_degrees(
        quotient.n,
        quotient.d,
        S,
        handoff_params,
        quotient,
    )
    @test all(
        degree -> degree <= 0 || degree in DeRham.cached_degrees(cache),
        expected_forward,
    )
    @test all(
        degree -> degree <= 0 || degree in DeRham.cached_reverse_lookups(cache),
        expected_reverse,
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
    empty_S_mixed = DeRham.reduction_setup(
        frontend,
        Int[],
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
    @test isempty(empty_S_mixed.S)
    @test empty_S_mixed.handoff_pole_order == n + 1
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
        dense_reduction = :akr,
    )
    @test_throws ArgumentError DeRham.reduction_setup(
        frontend,
        Int[],
        :depthfirst;
        dense_reduction = :akr,
    )
    @test_throws ArgumentError DeRham.reduction_setup(
        frontend,
        [-1],
        :depthfirst;
        dense_reduction = :akr_levelbylevel,
        handoff_pole_order = 1,
    )
end
