
function reducetransform(FT, N_m, S, f, P, pseudoInverseMat, p, params, cache, context)
    if params.algorithm == :pchunk
        reducetransform_pchunk(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    elseif params.algorithm == :depthfirst
        reducetransform_depthfirst(
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
        )
    elseif params.algorithm == :varbyvar
        reducetransform_varbyvar(FT, N_m, S, f, pseudoInverseMat, p, cache, params, context)
    elseif params.algorithm == :akr
        reducetransform_akr(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    else
        println(params.algorithm)
        throw(ArgumentError("Unsupported Algorithm: "))
    end
end

function reducetransform(FT, N_m, S, f, pseudoInverseMat, p, params, cache, context)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    P = HypersurfacePresentation(n, d)
    return reducetransform(FT, N_m, S, f, P, pseudoInverseMat, p, params, cache, context)
end
