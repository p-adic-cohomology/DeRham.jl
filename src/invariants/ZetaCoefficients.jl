"""
    zeta_coefficients_from_frobenius(FM, n, q, polygon, relative_precision, verbose)

Given the Frobenius matrix, computes the corresponding L-polynomial det(1-tq^{-1}FM)
returns an Oscar polynomial if as_poly=true, and the coefficients if as_poly=false

INPUTS:
* "FM" -- Frobenius matrix
* "n" -- integer, dimension of the ambient projective space
* "q" -- integer, characteristic of base field

"""

function zeta_coefficients_from_frobenius(FM, n, q, polygon, relative_precision, verbose)
    @assert size(FM, 1) == size(FM, 2) "FM is not a square matrix"

    P, T = polynomial_ring(ZZ, "T")
    lift_to_int(s) = map(x -> lift(ZZ, x), s)

    f = charpoly(P, lift_to_int(FM))
    cp_coeffs = collect(coefficients(f))
    return compute_Lpolynomial(n, q, polygon, relative_precision, cp_coeffs, verbose)
end


"""
    zeta_coefficients(f; S=[-1], verbose=0, changef=true, givefrobmat=false,
                      algorithm=:default, use_gpu=false, use_threads=false)

Computes zeta coefficients for the projective hypersurface defined by the
homogeneous Oscar polynomial `f` over a finite field of characteristic `p`.
The returned vector contains the coefficients of the middle-cohomology
L-polynomial in descending powers of `T`, as `ZZRingElem`s, with constant
coefficient `1` last. If `f` is not smooth, the polynomial-based method
returns `false`.

Keyword arguments:
- `S=[-1]`: controls the S-smooth model used by the reduction; the default
  sentinel expands to all indices `0:n`, where `n` is the dimension of the
  ambient projective space.
- `verbose=0`: controls diagnostic output; zero prints nothing and larger
  values print progressively more detail.
- `changef=true`: if the input is not an S-smooth model, allows the algorithm
  to try a random change of variables to find one. With `false`, or if no
  suitable model is found, the method returns `false`.
- `givefrobmat=false`: if true, returns `(F, coeffs)`, the Frobenius matrix
  and the coefficient vector, instead of just the coefficient vector.
- `algorithm=:default`: selects the controlled-reduction algorithm; `:default`
  chooses automatically.
- `use_gpu=false`: enables GPU execution for supported operations.
- `use_threads=false`: enables threaded execution for supported operations.

The `termorder`, `vars_reversed`, `fastevaluation`, `always_use_bigints`,
`context`, and `precision_info` keywords are internal or experimental and
may change. Their current defaults are `:invlex`, `false`, `true`, `false`,
`nothing`, and `nothing`, respectively.

The returned coefficients can be passed to [`L_polynomial`](@ref),
[`zeta_function`](@ref), and other coefficient-taking invariants.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> DeRham.zeta_coefficients(f)
3-element Vector{ZZRingElem}:
 7
 1
 1
```
"""
function zeta_coefficients(
    f;
    S = [-1],
    verbose = 0,
    changef = true,
    givefrobmat = false,
    algorithm = :default,
    termorder = :invlex,
    vars_reversed = false,
    fastevaluation = true,
    always_use_bigints = false,
    use_gpu = false,
    use_threads = false,
    context = nothing,
    precision_info = nothing,
)
    PR = parent(f)
    R = coefficient_ring(PR)
    p = Int64(characteristic(PR))
    q = p
    n = nvars(PR) - 1
    d = total_degree(f)

    if S == [-1]
        S = collect(0:n)
    end

    if length(S) > d
        error("Length of S must be <= $d")
    end

    if algorithm == :default
        if S == [n]
            algorithm = :varbyvar
        elseif n > 4
            for i = 0:n
                if is_Ssmooth(f, [n-i])
                    algorithm = :varbyvar
                    S = [n]
                    vars = gens(PR)
                    new_vars = copy(vars)
                    for j = 0:n
                        if j == i
                            new_vars[i+1] = vars[n+1]
                        end
                    end
                    new_vars[n+1] = vars[i+1]
                    break
                elseif i == n
                    algorithm = :depthfirst
                end
            end
        else
            algorithm = :depthfirst
        end
    end

    if algorithm==:varbyvar && (S != [n])
        throw("S must be [$n] for varbyvar")
    end

    # vars_reversed = false
    params = ZetaFunctionParams(
        verbose,
        givefrobmat,
        algorithm,
        termorder,
        vars_reversed,
        fastevaluation,
        always_use_bigints,
        use_gpu,
        use_threads,
    )

    cache = controlled_reduction_cache(n, d, S, params)

    (0 < verbose) && println("p = $p")
    (9 < verbose) && println("Working with a degree $d hypersurface in P^$n")

    (basis, Basis) = get_basis_of_cohomology_twoflavors(f, S, params, cache)

    #println("Basis of cohomology is $Basis")
    (9 < verbose) && println("Basis of cohomology is $Basis")

    if precision_info === nothing
        (hodge_polygon, r_m, N_m, M) = calculate_precision_information(f, Basis, verbose)
    else
        (hodge_polygon, r_m, N_m, M) = precision_info
    end

    (1 < verbose) && println("N_m=$N_m")

    (9 < verbose) && println(
        "We work modulo $p^$M, and compute up to the $N_m-th term of the Frobenius power series",
    )

    if always_use_bigints || BigInt(2)^64 < BigInt(p)^M
        residue = BigInt(p)^M
    elseif BigInt(2)^63 < BigInt(p)^M
        # If we use UInt, we get things between 2^63 and 2^64
        residue = UInt(p)^M
    else
        residue = p^M
    end

    # FOR DEBUGGING
    #residue = BigInt(p)^M

    precisionring, = residue_ring(ZZ, residue)
    precisionringpoly, pvars = polynomial_ring(precisionring, ["x$i" for i = 0:n])

    if (0 < verbose)
        println("Computing the relations matrix...")
        @time f_changed, f, pseudo_inverse_mat_new =
            find_Ssmooth_model(f, M, S, params, changef, cache)
    else
        f_changed, f, pseudo_inverse_mat_new =
            find_Ssmooth_model(f, M, S, params, changef, cache)
    end

    # For large examples, force a GC here; Julia's GC doesn't seem to automatically cause the
    # garbage collector to run if GPU allocations become too much.
    if 5 <= n
        if (0 < verbose)
            println("collecting garbage...")
            @time GC.gc()
        else
            GC.gc()
        end
    end

    if (f == false)
        (0 < verbose) && println("f is not smooth (or not S-smooth) and we're done. ")
        return false
    end

    (9 < verbose) && println("pseudo_inverse_mat is $pseudo_inverse_mat_new")

    # recomputes basis if f is different
    if f_changed
        (0 < verbose) && println("New model is $f")
        basis = compute_monomial_bases(f, params, cache) # basis of cohomology
        Basis = []
        for i = 1:n
            for j in basis[i]
                push!(Basis, [j, i])
            end
        end
        (9 < verbose) && println("New basis of cohomology is $Basis")
    end

    for i = 1:length(Basis)
        Basis[i][1] = liftCoefficients(precisionring, precisionringpoly, Basis[i][1])
    end

    (9 < verbose) && println("Basis of cohomology is $Basis")


    fLift = liftCoefficients(precisionring, precisionringpoly, f)
    if (0 < verbose)
        println("\nApplying Frobenius to basis elements")
        @time FBasis = applyFrobeniusToBasis(Basis, fLift, N_m, p, params)
    else
        FBasis = applyFrobeniusToBasis(Basis, fLift, N_m, p, params)
    end

    l = d * n - n + d - length(S)


    MS = matrix_space(
        precisionring,
        nrows(pseudo_inverse_mat_new),
        ncols(pseudo_inverse_mat_new),
    )
    pseudo_inverse_mat = MS()

    if (0 < verbose)
        println("Computing T matrix...")
        @time T = computeT(f, basis, M, params, cache)
    else
        T = computeT(f, basis, M, params, cache)
    end
    (9 < verbose) && println("T matrix is $T")

    for i = 1:nrows(pseudo_inverse_mat_new)
        for j = 1:ncols(pseudo_inverse_mat_new)
            pseudo_inverse_mat[i, j] = ZZ(pseudo_inverse_mat_new[i, j])
        end
    end

    #TODO: check which algorithm we're using
    (0 < verbose) && println("\nStarting controlled reduction...")
    Reductions = reducetransform(
        FBasis,
        N_m,
        S,
        fLift,
        pseudo_inverse_mat,
        p,
        params,
        cache,
        context,
    )
    (2 < verbose) && println(Reductions)

    ev = cache[n*d-n-1]
    (9 < verbose) && println(convert_p_to_m([Reductions[1][1][1], Reductions[2][1][1]], ev))
    FM = compute_frobenius_matrix(n, p, d, N_m, Reductions, T, Basis, params, cache)
    (9 < verbose) && println("The Frobenius matrix is $FM")

    if givefrobmat
        (FM, zeta_coefficients_from_frobenius(FM, n, q, hodge_polygon, r_m, verbose))
    else
        zeta_coefficients_from_frobenius(FM, n, q, hodge_polygon, r_m, verbose)
    end
end

"""
a wrapper to zeta_coefficients

INPUTS:
* "f" -- Oscar polynomial (should be homogeneous) over the integers
* "p" -- a prime number, integer

"""
function zeta_coefficients(
    f,
    p;
    S = [-1],
    verbose = 0,
    changef = true,
    givefrobmat = false,
    algorithm = :default,
    termorder = :invlex,
    vars_reversed = false,
    fastevaluation = false,
    always_use_bigints = false,
    use_gpu = false,
    use_threads = false,
    context = nothing,
)
    @assert is_prime(p) "p must be prime"
    PR = parent(f)
    PRmodp, hom = change_base_ring(GF(p), PR)

    if algorithm == :default
        if S == [n]
            algorithm = :varbyvar
        elseif n > 4
            for i = 0:n
                if is_Ssmooth(f, [n-i])
                    algorithm = :varbyvar
                    S = [n]
                    vars = gens(PR)
                    new_vars = copy(vars)
                    for j = 0:n
                        if j == i
                            new_vars[i+1] = vars[n+1]
                        end
                    end
                    new_vars[n+1] = vars[i+1]
                    break
                elseif i == n
                    algorithm = :depthfirst
                end
            end
        else
            algorithm = :depthfirst
        end
    end

    return zeta_coefficients(
        hom(f);
        S = S,
        verbose = verbose,
        changef = changef,
        givefrobmat = givefrobmat,
        algorithm = algorithm,
        termorder = termorder,
        vars_reversed = vars_reversed,
        fastevaluation = fastevaluation,
        always_use_bigints = always_use_bigints,
        use_gpu = use_gpu,
        use_threads = use_threads,
        context = context,
    )
end


"""
    zeta_coefficients_with_precision(f, r::Integer; basis=nothing, verbose=0, kwargs...)

Compute the zeta coefficients for the hypersurface defined by `f`, using `r`
as a requested relative `p`-adic precision for every Hodge component. This
scalar overload applies the same request to all components. The algorithm
combines the requested precision with its automatically computed minimum.

- `f`: homogeneous polynomial over a finite field.
- `r`: requested relative precision, repeated across the Hodge components.
- `basis`: optional cohomology basis; the package computes one by default.
- `verbose`: diagnostic verbosity.
- `kwargs...`: keywords forwarded to [`zeta_coefficients`](@ref).

Returns the descending `ZZRingElem` coefficient vector, or `false` when no
smooth model is found.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> length(DeRham.zeta_coefficients_with_precision(f, 1))
3
```
"""
function zeta_coefficients_with_precision(
    f,
    r::Integer;
    basis = nothing,
    verbose = 0,
    kwargs...,
)
    r_m = fill(Int(r), nvars(parent(f)) - 1)
    return zeta_coefficients_with_precision(
        f,
        r_m;
        basis = basis,
        verbose = verbose,
        kwargs...,
    )
end

"""
    zeta_coefficients_with_precision(f, r_m::AbstractVector; basis=nothing, verbose=0, kwargs...)

Vector-precision overload of [`zeta_coefficients_with_precision`](@ref).
`r_m` gives one requested relative precision for each Hodge component and
must have the length required by the basis. The algorithm takes the
componentwise maximum of these requests and its automatically computed
minimum. The return value and other keywords are the same as for the scalar
overload.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> length(DeRham.zeta_coefficients_with_precision(f, [1, 1]))
3
```
"""
function zeta_coefficients_with_precision(
    f,
    r::AbstractVector;
    basis = nothing,
    verbose = 0,
    kwargs...,
)
    precision_info =
        precision_information_max_auto_r_m(f, r; basis = basis, verbose = verbose)
    return zeta_coefficients(
        f;
        precision_info = precision_info,
        verbose = verbose,
        kwargs...,
    )
end

"""
    zeta_coefficients_with_precision(f, hodge_numbers, r_m, N_m, M; basis=nothing, verbose=0, kwargs...)

Advanced overload that supplies the Hodge numbers and the algorithm's
precision data directly. These parameters are tied to the current algorithm
and may change; prefer one of the `r` overloads unless you need to reproduce
the current precision calculation. Returns the same value as
[`zeta_coefficients`](@ref).
"""
function zeta_coefficients_with_precision(
    f,
    hodge_numbers::AbstractVector,
    r_m::AbstractVector,
    N_m::AbstractVector,
    M::Integer;
    basis = nothing,
    verbose = 0,
    kwargs...,
)
    precision_info = (SlopesPolygon(hodge_numbers), r_m, N_m, M)
    return zeta_coefficients(
        f;
        precision_info = precision_info,
        verbose = verbose,
        kwargs...,
    )
end
