# TODO: PR this into Oscar, hopefully

function artin_shreier_curve(num, denom)
    p = characteristic(parent(num))
    R, (x, y) = polynomial_ring(GF(p), ["x", "y"])

    lhs = evaluate(num, x)
    prod_factor = evaluate(denom, x)

    #TODO: make this into a scheme in Oscar's language
    prod_factor * (y^p - y) - lhs
end

function projective_artin_shreier_curve(num, denom)
    nonhomog = artin_shreier_curve(num, denom)
    P = parent(nonhomog)
    H = homogenizer(P, "z")
    H(nonhomog)
end

#function artih_shreier_curve(p,f)
#  TODO:
#
#end

function artin_shreier_witt_curve(p, num, denom)
    #TODO: make heuristic algorithm for delta_1 in MMPSingularities.jl

    error("Artin-Shreier-Witt curves not implemented yet")
end

"""
    fermat_hypersurface(n, d, p)

Construct the Fermat hypersurface polynomial ``x_1^d + \\cdots + x_n^d`` over the
finite field ``\\mathbb{F}_p``. It has ``n`` variables and defines a hypersurface of
dimension ``n - 2`` in projective space.

- `n`: number of variables.
- `d`: degree.
- `p`: prime characteristic of the coefficient field.

Returns the homogeneous Oscar polynomial.

### Examples

```jldoctest
julia> f = DeRham.fermat_hypersurface(3, 3, 7);

julia> total_degree(f)
3
```
"""
function fermat_hypersurface(n, d, p)
    F = GF(p)
    R, _ = polynomial_ring(F, n)


    C = MPolyBuildCtx(R)
    exp_vec = zeros(Int, n)
    for i = 1:n
        zero!(exp_vec)
        exp_vec[i] = d
        push_term!(C, F(1), exp_vec)
    end

    finish(C)
end

"""
Gives the hypersurface of degree d and dimension n-2
which is defined by the polynomial which has every
term with coefficient 1.
"""
function full_hypersurface(n, d, p)
    F = GF(p)
    R, _ = polynomial_ring(F, n)

    evs = gen_exp_vec(n, d)

    C = MPolyBuildCtx(R)

    for ev in evs
        push_term!(C, F(1), ev)
    end

    finish(C)
end

"""
    random_hypersurface(n, d, p)

Construct a homogeneous degree-``d`` polynomial in ``n`` variables over ``\\mathbb{F}_p``,
with independent uniform random coefficients for its degree-``d`` terms. It
defines a hypersurface of dimension ``n - 2`` in projective space.

- `n`: number of variables.
- `d`: degree of the terms.
- `p`: prime characteristic of the coefficient field.

Returns the Oscar polynomial. The result is random and is not guaranteed to
define a smooth hypersurface.

### Examples

```jldoctest
julia> f = DeRham.random_hypersurface(3, 3, 7);

julia> nvars(parent(f))
3
```
"""
function random_hypersurface(n, d, p)
    F = GF(p)
    R, _ = polynomial_ring(F, n)

    evs = gen_exp_vec(n, d)
    cs = rand(0:(p-1), length(evs))

    C = MPolyBuildCtx(R)

    for i in eachindex(evs)
        push_term!(C, F(cs[i]), evs[i])
    end

    finish(C)
end
