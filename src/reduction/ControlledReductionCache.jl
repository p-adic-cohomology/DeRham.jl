_reduction_g_degree(P::CohomologyQuotient) = numerator_degree(P, terminal_pole_order(P)) + 1

function _reduction_handoff_order(P::CohomologyQuotient, params)
    order = something(params.handoff_pole_order, terminal_pole_order(P))
    order >= terminal_pole_order(P) ||
        throw(ArgumentError("handoff pole order is below the terminal pole order"))
    return order
end

_reduction_g_degree(P::CohomologyQuotient, params) =
    numerator_degree(P, _reduction_handoff_order(P, params)) + 1

"""
    controlled_reduction_cache

Give the minimal PolyExpCache needed for controlled reduction
as in Prop 1.15 of Costa's thesis.

n in this function is the number of variables minus one.

For each chain degree (the original terminal degree and, when different, the
dense handoff degree), with `s = length(S)`, the degrees needed are

FORWARD:
For h in pole_orders(P), numerator_degree(P, h) and numerator_degree(P, h) - (d - 1), when nonnegative
g
g - 1
g - s
g - s + 1
g + d - s
d

REVERSE:
d
g
g + d - s

fields
------
n - int
d - int
S - vector of ints
params - the ControlledReductionParamaters
P - cohomology quotient; defaults to the hypersurface quotient for n and d
"""
function controlled_reduction_cache(n, d, S, params, P = HypersurfaceQuotient(n, d))
    degsforward, degsreverse = _controlled_reduction_cache_degrees(n, d, S, params, P)
    PolyExpCache(n+1, params.termorder, degsforward, degsreverse, vars_reversed = false)
end

function _controlled_reduction_cache_degrees(n, d, S, params, P)
    chain_degrees = unique([_reduction_g_degree(P), _reduction_g_degree(P, params)])
    s = length(S)
    degsforward = Int[]
    degsreverse = Int[]
    for g in chain_degrees
        append!(degsforward, [g, g - 1, g - s, g - s + 1, g + d - s, d])
        append!(degsreverse, [d, g, g + d - s])
    end
    for h in pole_orders(P)
        numerator = numerator_degree(P, h)
        if 0 ≤ numerator
            push!(degsforward, numerator)
        end
        denominator = numerator - (d - 1)
        if 0 ≤ denominator
            push!(degsforward, denominator)
        end
    end

    return degsforward, degsreverse
end

function extend_controlled_reduction_cache!(cache::PolyExpCache, n, d, S, params, P)
    cache.n == n + 1 || throw(ArgumentError("cache has the wrong number of variables"))
    cache.termorder == params.termorder || throw(ArgumentError("cache term order differs"))
    cache.vars_reversed == false || throw(ArgumentError("cache variable order differs"))

    degsforward, degsreverse = _controlled_reduction_cache_degrees(n, d, S, params, P)
    for degree in degsforward
        generate_degree_forward(cache, degree)
    end
    for degree in degsreverse
        generate_degree_reverse(cache, degree)
    end
    return cache
end
