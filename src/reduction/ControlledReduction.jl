
function reducetransform(
    FT,
    N_m,
    S,
    f,
    P,
    fe::ZetaFrontend,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    setup === nothing && (
        setup = reduction_setup(
            fe,
            S,
            params.algorithm;
            use_gpu = params.use_gpu,
            use_threads = params.use_threads,
            always_use_bigints = params.always_use_bigints,
            fastevaluation = params.fastevaluation,
            context = context,
        )
    )
    if setup.policy == :pchunk
        reducetransform_pchunk(FT, N_m, S, f, P, pseudoInverseMat, p, cache, params, fe)
    elseif setup.policy == :depthfirst
        reducetransform_depthfirst(
            FT,
            N_m,
            S,
            f,
            P,
            fe,
            pseudoInverseMat,
            p,
            cache,
            params,
            context,
            setup,
        )
    elseif setup.policy == :varbyvar
        reducetransform_varbyvar(
            FT,
            N_m,
            S,
            f,
            P,
            pseudoInverseMat,
            p,
            cache,
            params,
            context,
            fe,
            setup,
        )
    elseif setup.policy == :akr
        reducetransform_akr(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    else
        println(setup.policy)
        throw(ArgumentError("Unsupported Algorithm: "))
    end
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    P::HypersurfaceQuotient,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        ProjectiveHypersurface(f),
        pseudoInverseMat,
        p,
        params,
        cache,
        context;
        setup = setup,
    )
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    P::DoubleCoverQuotient,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        DoubleCover(f),
        pseudoInverseMat,
        p,
        params,
        cache,
        context;
        setup = setup,
    )
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    pseudoInverseMat,
    p,
    params,
    cache,
    context;
    setup = nothing,
)
    fe = ProjectiveHypersurface(f)
    P = cohomology_quotient(fe)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        fe,
        pseudoInverseMat,
        p,
        params,
        cache,
        context;
        setup = setup,
    )
end

function reducetransform(FT, N_m, S, f, pseudoInverseMat, p, params, cache; setup = nothing)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        pseudoInverseMat,
        p,
        params,
        cache,
        nothing;
        setup = setup,
    )
end
