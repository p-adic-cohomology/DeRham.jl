@testset "reduction selection characterization" begin
    R, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    cubic = x^3 + y^3 + z^3
    P = DeRham.HypersurfaceQuotient(2, 3)

    for (S, expected_algorithm) in
        (([0], :depthfirst), ([2], :varbyvar), ([0, 1, 2], :depthfirst))
        algorithm, selected_S, selected_variable =
            DeRham._default_reduction_selection(cubic, P, S)
        @test algorithm == expected_algorithm
        @test selected_S == S
        @test selected_variable === nothing
    end

    R5, variables = polynomial_ring(GF(7), ["x$i" for i = 0:5])
    quadric = sum(variable^2 for variable in variables)
    P5 = DeRham.HypersurfaceQuotient(5, 2)
    algorithm, selected_S, selected_variable =
        DeRham._default_reduction_selection(quadric, P5, [0, 1])
    @test algorithm == :varbyvar
    @test selected_S == [5]
    @test selected_variable == 0

    singular_quadric = variables[1]^2
    algorithm, selected_S, selected_variable =
        DeRham._default_reduction_selection(singular_quadric, P5, [0, 1])
    @test algorithm == :depthfirst
    @test selected_S == [0, 1]
    @test selected_variable === nothing

    automatic_input = sum(variables[i]^2 for i = 1:5) + variables[4] * variables[6]
    automatic_selection_error = try
        DeRham.zeta_coefficients(
            automatic_input;
            S = [0],
            algorithm = :default,
            changef = false,
            precision_info = (nothing, [1], [1], 4),
        )
        nothing
    catch err
        err
    end
    @test automatic_selection_error isa ArgumentError
    if automatic_selection_error isa ArgumentError
        @test occursin(
            "matrix from f is not right invertible",
            sprint(showerror, automatic_selection_error),
        )
    end

    hypersurface = DeRham.ProjectiveHypersurface(cubic)
    unsupported_algorithm_error = try
        DeRham.zeta_coefficients(hypersurface; algorithm = :unsupported)
        nothing
    catch err
        err
    end
    @test unsupported_algorithm_error isa ArgumentError
    if unsupported_algorithm_error isa ArgumentError
        @test occursin(
            "not yet supported for",
            sprint(showerror, unsupported_algorithm_error),
        )
    end
    @test DeRham.varbyvar_reduction_setup(hypersurface, [2]) ==
          (:varbyvar, [2], nothing, [2])
    @test DeRham.varbyvar_reduction_setup(hypersurface, [0]) ==
          (:depthfirst, [0], nothing, [2])

    Rdc, (s, t) = polynomial_ring(GF(5), ["s", "t"])
    double_cover = DeRham.DoubleCover(s^4 + t^4)
    double_cover_terminal = DeRham.cohomology_quotient(double_cover).n
    @test DeRham._supports_reduction_algorithm(hypersurface, :varbyvar)
    @test DeRham._supports_reduction_algorithm(double_cover, :varbyvar)
    @test !DeRham._supports_reduction_algorithm(double_cover, :pchunk)
    @test DeRham.varbyvar_reduction_setup(double_cover, [double_cover_terminal]) ==
          (:depthfirst, [double_cover_terminal], nothing, [double_cover_terminal])

    setup_modulus = BigInt(5)^8
    setup_ring, = residue_ring(ZZ, setup_modulus)
    _, setup_variables = polynomial_ring(setup_ring, ["x0", "x1", "x2"])
    setup_f = sum(variable^3 for variable in setup_variables)
    setup_frontend = DeRham.ProjectiveHypersurface(setup_f)
    setup_P = DeRham.cohomology_quotient(setup_frontend)
    setup_configuration =
        DeRham.reduction_setup(setup_frontend, [2], :depthfirst; use_threads = true)
    setup_params = DeRham.ZetaFunctionParams(
        0,
        false,
        :depthfirst,
        :invlex,
        false,
        true,
        false,
        false,
        true,
    )
    setup_size = binomial(DeRham._reduction_g_degree(setup_P) + setup_P.n, setup_P.n)
    setup_matspace = matrix_space(setup_ring, setup_size, setup_size)
    setup_cache = DeRham.controlled_reduction_cache(2, 3, [2], setup_params, setup_P)
    execution_setup = DeRham.reduction_execution_setup(
        setup_configuration,
        setup_frontend,
        setup_f,
        setup_P,
        setup_matspace(),
        5,
        setup_cache,
        setup_params,
        nothing,
    )
    @test setup_configuration.allowed_policies == (:depthfirst, :pchunk, :varbyvar, :akr)
    @test setup_configuration.lazy_Ruv
    @test execution_setup.pep_type == DeRham.LazyPEP{typeof(setup_matspace())}
    @test execution_setup.backend == :cpu_uint
    @test execution_setup.representation == :uint
    @test execution_setup.kernel === DeRham.eval_to_linear!
    @test execution_setup.threads
    mismatched_thread_params = DeRham.ZetaFunctionParams(
        setup_params.verbose,
        setup_params.givefrobmat,
        setup_params.algorithm,
        setup_params.termorder,
        setup_params.vars_reversed,
        setup_params.fastevaluation,
        setup_params.always_use_bigints,
        setup_params.use_gpu,
        false,
    )
    @test_throws ArgumentError DeRham._validate_reduction_execution_setup(
        setup_configuration,
        setup_frontend,
        mismatched_thread_params,
    )
    rebound_frontend = DeRham.ProjectiveHypersurface(setup_f + setup_variables[1]^3)
    @test DeRham._validate_reduction_execution_setup(
        setup_configuration,
        rebound_frontend,
        setup_params,
    ) === nothing
    mismatched_frontend = DeRham.ProjectiveHypersurface(setup_f * setup_f)
    @test_throws ArgumentError DeRham._validate_reduction_execution_setup(
        setup_configuration,
        mismatched_frontend,
        setup_params,
    )

    cases = ((2, 3, [0]), (2, 3, [2]), (2, 3, [0, 1, 2]), (4, 3, [0]), (4, 4, [0]))
    moduli = (BigInt(5)^8, BigInt(2)^40, BigInt(2)^110)

    for (n, d, S) in cases
        _, variables = polynomial_ring(GF(5), ["x$i" for i = 0:n])
        f = sum(variable^d for variable in variables)
        P = DeRham.HypersurfaceQuotient(n, d)
        expected_algorithm = S == [n] ? :varbyvar : :depthfirst
        algorithm, selected_S, _ = DeRham._default_reduction_selection(f, P, S)
        @test algorithm == expected_algorithm
        @test selected_S == S

        for modulus in moduli,
            use_gpu in (false, true),
            use_threads in (false, true),
            always_use_bigints in (false, true),
            fastevaluation in (false, true)

            params = DeRham.ZetaFunctionParams(
                0,
                false,
                algorithm,
                :invlex,
                false,
                fastevaluation,
                always_use_bigints,
                use_gpu,
                use_threads,
            )
            ring, = residue_ring(ZZ, modulus)
            matspace = matrix_space(ring, 1, 1)
            cache = DeRham.controlled_reduction_cache(n, d, S, params, P)
            lazy = use_gpu || length(S) < d || d < n || algorithm == :varbyvar
            compute = _ -> [matspace() for _ = 1:(n+2)]
            pep = DeRham.select_Ruv_PEP(n, d, S, params, compute, lazy, matspace, cache)

            if use_gpu
                if BigInt(2)^25 < modulus < BigInt(2)^106
                    @test typeof(pep) == DeRham.LazyPEP{DeRham.KaratsubaMatrix{Float64}}
                else
                    @test typeof(pep) == DeRham.LazyPEP{DeRham.CuModMatrix{Float64}}
                end
                @test DeRham.default_context_type(matspace, params) ==
                      DeRham.ControlledReductionContext{
                    DeRham.CuModMatrix{Float64},
                    DeRham.CuModVector{Float64},
                }
            else
                if lazy
                    @test typeof(pep) == DeRham.LazyPEP{typeof(matspace())}
                else
                    @test typeof(pep) == DeRham.EagerPEP{typeof(matspace())}
                end
                context = DeRham.default_context(matspace, pep, params)
                @test typeof(context.A) == typeof(matspace())
                if always_use_bigints || modulus > BigInt(2)^64
                    @test DeRham.default_context_type(matspace, params) ==
                          DeRham.ControlledReductionContext{
                        DeRham.ZZModMatrix,
                        Vector{DeRham.ZZRingElem},
                    }
                    @test typeof(context.g) == Vector{DeRham.ZZRingElem}
                else
                    @test DeRham.default_context_type(matspace, params) ==
                          DeRham.ControlledReductionContext{DeRham.zzModMatrix,Vector{UInt}}
                    @test typeof(context.g) == Vector{UInt}
                end
            end

            kernel = DeRham._varbyvar_evaluation_kernel(params, d, n, modulus)
            if use_gpu && !(BigInt(2)^25 < modulus < BigInt(2)^106)
                @test kernel === DeRham.eval_to_linear_gpu!
            elseif use_gpu && d == 3 && (n == 4 || n == 5)
                @test kernel === DeRham.eval_to_linear_gpu_karatsuba!
            else
                @test kernel === DeRham.eval_to_linear!
            end
        end
    end
end
