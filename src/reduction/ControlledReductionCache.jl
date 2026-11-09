_reduction_g_degree(P::CohomologyPresentation) =
    numerator_degree(P, terminal_pole_order(P)) + 1

"""
    controlled_reduction_cache

Give the minimal PolyExpCache needed for controlled reduction
as in Prop 1.15 of Costa's thesis.

n in this function is the number of variables minus one.

With `g = _reduction_g_degree(P)` and `s = length(S)`, the degrees needed are

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
P - cohomology presentation; defaults to the hypersurface presentation for n and d
"""
function controlled_reduction_cache(n, d, S, params, P = HypersurfacePresentation(n, d))
    g = _reduction_g_degree(P)
    s = length(S)
    degsforward = [
        g,
        g - 1,
        g - s,
        g - s + 1,
        g + d - s,
        d,
    ]
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

    degsreverse = [d, g, g + d - s]

    PolyExpCache(n+1, params.termorder, degsforward, degsreverse, vars_reversed = false)

end
