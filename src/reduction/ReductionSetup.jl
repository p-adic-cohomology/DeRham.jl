struct ReductionSetup{F,S}
    frontend::F
    policy::Symbol
    S::S
    selected_variable::Union{Nothing,Int}
    required_S::Vector{Int}
    allowed_policies::Tuple
    lazy_Ruv::Bool
    use_gpu::Bool
    use_threads::Bool
    always_use_bigints::Bool
    fastevaluation::Bool
end

struct ReductionExecutionSetup{C}
    configuration::ReductionSetup
    Ruv::Any
    ruv_function::Any
    context_factory::C
    context_type::Any
    kernel::Function
    pep_type::Any
    backend::Symbol
    representation::Symbol
    threads::Bool
    context_provided::Bool
    preconfigured_engine::Any
end

function _validate_reduction_setup(fe, algorithm; use_gpu = false, context = nothing)
    if !_supports_reduction_algorithm(fe, algorithm)
        throw(ArgumentError("algorithm $algorithm is not yet supported for $(typeof(fe))"))
    end

    if !(fe isa ProjectiveHypersurface)
        if use_gpu
            throw(ArgumentError("GPU execution is not yet supported for $(typeof(fe))"))
        elseif context !== nothing
            throw(
                ArgumentError(
                    "custom reduction contexts are not yet supported for $(typeof(fe))",
                ),
            )
        end
    end
    return nothing
end

"""
    reduction_setup(fe, S, algorithm; kwargs...)

Resolve the frontend's allowed policy, automatic policy and variable selection,
and policy-specific R_uv storage strategy. `reduction_execution_setup` completes
the runtime choices once the precision ring is available. The runtime setup
preserves the existing PEP and evaluation-kernel selectors, including their
different GPU Karatsuba eligibility rules.
"""
function reduction_setup(
    fe::ZetaFrontend,
    S,
    algorithm;
    use_gpu = false,
    use_threads = false,
    always_use_bigints = false,
    fastevaluation = true,
    context = nothing,
)
    allowed_policies = _reduction_allowed_policies(fe)
    _validate_reduction_setup(fe, algorithm; use_gpu = use_gpu, context = context)

    if algorithm == :default
        policy, selected_S, selected_variable, required_S = varbyvar_reduction_setup(fe, S)
    else
        policy = algorithm
        selected_S = S
        selected_variable = nothing
        required_S = [cohomology_quotient(fe).n]
    end

    if policy == :varbyvar && selected_S != required_S
        throw("S must be $required_S for varbyvar")
    end

    P = cohomology_quotient(fe)
    d = P.d
    n = P.n
    lazy_Ruv = if policy == :varbyvar
        true
    elseif policy == :depthfirst
        length(selected_S) < d || d < n
    else
        policy == :pchunk
    end

    return ReductionSetup(
        fe,
        policy,
        selected_S,
        selected_variable,
        required_S,
        allowed_policies,
        lazy_Ruv,
        use_gpu,
        use_threads,
        always_use_bigints,
        fastevaluation,
    )
end

function _reduction_backend(params, m)
    if params.use_gpu
        return ZZ(2)^25 < m < ZZ(2)^106 ? :gpu_karatsuba : :gpu_cumod
    elseif params.always_use_bigints || ZZ(2)^64 < m
        return :cpu_bigint
    end
    return :cpu_uint
end

function _reduction_representation(vector)
    if vector isa Vector{ZZRingElem}
        return :bigint
    elseif vector isa Vector{UInt}
        return :uint
    elseif vector isa CuModVector
        return :cumod
    elseif vector isa KaratsubaVector
        return :karatsuba
    end
    return :custom
end

function _validate_reduction_execution_setup(configuration, fe, params)
    frontend_matches =
        typeof(fe) === typeof(configuration.frontend) &&
        cohomology_quotient(fe) == cohomology_quotient(configuration.frontend)
    frontend_matches ||
        throw(ArgumentError("the frontend does not match the reduction setup"))
    params.algorithm == configuration.policy ||
        throw(ArgumentError("the algorithm does not match the reduction setup"))
    params.use_gpu == configuration.use_gpu ||
        throw(ArgumentError("the GPU option does not match the reduction setup"))
    params.use_threads == configuration.use_threads ||
        throw(ArgumentError("the threading option does not match the reduction setup"))
    params.always_use_bigints == configuration.always_use_bigints || throw(
        ArgumentError(
            "the integer representation option does not match the reduction setup",
        ),
    )
    params.fastevaluation == configuration.fastevaluation ||
        throw(ArgumentError("the evaluation option does not match the reduction setup"))
    return nothing
end

function reduction_execution_setup(
    configuration::ReductionSetup,
    fe::ZetaFrontend,
    f,
    P,
    pseudoInverseMat,
    p,
    cache,
    params,
    context,
)
    _validate_reduction_execution_setup(configuration, fe, params)
    d = total_degree(f)
    n = nvars(parent(f)) - 1
    g = _reduction_g_degree(P)
    g_length = binomial(g + n, n)
    MS1 = matrix_space(coefficient_ring(parent(f)), g_length, g_length)
    m = Integer(modulus(base_ring(MS1)))

    ruv_function =
        reduction_ruv_function(fe, configuration.S, f, pseudoInverseMat, cache, params)
    computeRuv = if 3 < params.verbose
        V -> begin
            println("Computing Ruv for V = $V for the first time.")
            @time ruv_function(V)
        end
    else
        ruv_function
    end

    context_provided = context !== nothing
    if !context_provided
        if 0 < params.verbose
            println("Creating the Ruv PEP object...")
            @time Ruv = select_Ruv_PEP(
                n,
                d,
                configuration.S,
                params,
                computeRuv,
                configuration.lazy_Ruv,
                MS1,
                cache,
            )
        else
            Ruv = select_Ruv_PEP(
                n,
                d,
                configuration.S,
                params,
                computeRuv,
                configuration.lazy_Ruv,
                MS1,
                cache,
            )
        end
        context_factory = () -> default_context(MS1, Ruv, params)
        context_type = default_context_type(MS1, params)
        preconfigured_engine = nothing
    else
        if ZZ(2)^25 < m < ZZ(2)^106
            compute_gpu = V -> KaratsubaMat.(computeRuv(V))
        else
            compute_gpu = V -> cuMod.(computeRuv(V))
        end
        compute_float = V -> float_entries.(computeRuv(V))
        if 3 < n && d == 3 && configuration.S == [n] && params.use_gpu
            bind_ruv_function!(context.Ruvs, compute_float)
        elseif params.use_gpu
            bind_ruv_function!(context.Ruvs, compute_gpu)
        else
            bind_ruv_function!(context.Ruvs, computeRuv)
        end
        Ruv = context.Ruvs
        preconfigured_engine = PreconfiguredReductionEngine(context)
        context_factory = () -> preconfigured_engine.context
        context_type = typeof(context)
    end

    kernel = if configuration.policy == :varbyvar
        _varbyvar_evaluation_kernel(params, d, n, m)
    else
        eval_to_linear!
    end
    backend = _reduction_backend(params, m)
    representation =
        context_provided ? _reduction_representation(context.g) :
        backend == :gpu_karatsuba ? :karatsuba :
        backend == :gpu_cumod ? :cumod : backend == :cpu_bigint ? :bigint : :uint

    return ReductionExecutionSetup(
        configuration,
        Ruv,
        ruv_function,
        context_factory,
        context_type,
        kernel,
        typeof(Ruv),
        backend,
        representation,
        !context_provided && configuration.use_threads,
        context_provided,
        preconfigured_engine,
    )
end
