
function is_Ssmooth(f, P, S)
    p = characteristic(parent(f))
    n = length(gens(parent(f))) - 1
    d = total_degree(f)
    l = _reduction_g_degree(P) + d - length(S)

    PR = parent(f)
    R = base_ring(PR)

    PRZZ, VarsZZ = polynomial_ring(ZZ, ["x$i" for i = 0:n])
    fLift = liftCoefficients(ZZ, PRZZ, f)

    params = default_params()
    cache = controlled_reduction_cache(n, d, S, params, P)

    M_ZZ = compute_controlled_matrix(fLift, l, S, ZZ, PRZZ, params, cache)

    M = matrix(R, M_ZZ)
    flag, B = Oscar.is_invertible_with_inverse(M, side = :right)
    flag
end

function is_Ssmooth(f, S)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    P = HypersurfaceQuotient(n, d)
    return is_Ssmooth(f, P, S)
end

function issmooth_linalg(f)
    n = length(gens(parent(f))) - 1
    is_Ssmooth(f, collect(0:n))
end

function random_change_of_variables(f)

    vars = gens(parent(f))
    n = length(vars)-1 # we don't need it here, but use the standard convention
    p = characteristic(parent(f))

    SLn = special_linear_group(n + 1, GF(p))
    mat = rand(SLn)
    new_vars = matrix(mat) * vars
    f_transformed = evaluate(f, new_vars)

    f_transformed
end

"""
    find_Ssmooth_model(f, M, S_target, params)
Find a change of variable that makes f S_target-smooth

Inputs:
* "f" -- Oscar polynomial (should be homogeneous)
* "M" -- integer, algorithm precision
* "S_target" -- list, target S-smoothness
* "params" -- the ControlledReductionParamaters
* "changef" -- bool
* "cache" -- the GradedExpCache used for this controlled reduction
"""
function find_Ssmooth_model(f, P, M, S_target, params, changef, cache)
    l = _reduction_g_degree(P) + total_degree(f) - length(S_target)
    return _find_Ssmooth_model(f, params, changef) do model
        pseudo_inverse_controlled_lifted(model, S_target, l, M, params, cache)
    end
end

function _find_Ssmooth_model(build_setup, f, params, changef)
    p = Int64(characteristic(parent(f)))
    PR = parent(f)
    n = nvars(PR) - 1
    SLn = special_linear_group(n + 1, GF(p))
    vars = gens(PR)

    f_transformed = f
    num_iter = min(2*p, 100)
    f_changed = false
    if !changef
        return f_changed, f_transformed, build_setup(f_transformed)
    end

    for i = 1:num_iter
        try
            return f_changed, f_transformed, build_setup(f_transformed)
        catch e
            (0 < params.verbose) &&
                println("This f is not S-smooth, changing to one that is")
            if isa(e, ArgumentError) && e.msg == "f is not smooth"

                return false, false, false
            else
                mat = rand(SLn)
                new_vars = matrix(mat) * vars
                f_transformed = evaluate(f, new_vars)
                f_changed = true
            end
        end
    end
    return (false, false, false)

end

function find_Ssmooth_model(f, M, S_target, params, changef, cache)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    P = HypersurfaceQuotient(n, d)
    return find_Ssmooth_model(f, P, M, S_target, params, changef, cache)
end
