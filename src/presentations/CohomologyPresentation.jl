"""
    abstract type CohomologyPresentation

A Griffiths–Dwork presentation of cohomology indexed by pole order.
"""
abstract type CohomologyPresentation end

"""
    HypersurfacePresentation(n, d)

The cohomology presentation for a degree-`d` hypersurface of dimension `n`.
"""
struct HypersurfacePresentation <: CohomologyPresentation
    n::Int
    d::Int
end

"""
    pole_orders(P::HypersurfacePresentation)

Return the pole orders that carry basis elements for the hypersurface `P`.
"""
pole_orders(P::HypersurfacePresentation) = 1:P.n

"""
    terminal_pole_order(P::HypersurfacePresentation)

Return the pole order where reduction of the hypersurface `P` stops.
"""
terminal_pole_order(P::HypersurfacePresentation) = P.n

"""
    numerator_degree(P::HypersurfacePresentation, m)

Return the degree of the numerator at pole order `m` for `P`.
"""
numerator_degree(P::HypersurfacePresentation, m::Integer) = P.d * m - (P.n + 1)

"""
    pole_exponent(P::HypersurfacePresentation, m)

Return the exponent of the hypersurface polynomial at pole order `m`.
"""
pole_exponent(::HypersurfacePresentation, m::Integer) = m

"""
    reduction_scalar(P::HypersurfacePresentation, m)

Return the exact Griffiths–Dwork scalar for the step from pole order `m + 1` to `m`.
"""
reduction_scalar(::HypersurfacePresentation, m::Integer) = QQ(1, m)

"""
    DoubleCoverPresentation(n, d)

The cohomology presentation for the double cover of `P^n` branched along a
degree-`d` hypersurface, with `d` even.
"""
struct DoubleCoverPresentation <: CohomologyPresentation
    n::Int
    d::Int

    function DoubleCoverPresentation(n::Integer, d::Integer)
        iseven(d) || throw(ArgumentError("the branch degree must be even"))
        new(Int(n), Int(d))
    end
end

"""
    pole_orders(P::DoubleCoverPresentation)

Return the pole orders that carry basis elements for the double cover `P`.
"""
pole_orders(P::DoubleCoverPresentation) = 1:(P.n + 1)

"""
    terminal_pole_order(P::DoubleCoverPresentation)

Return the pole order where reduction of the double cover `P` stops.
"""
terminal_pole_order(P::DoubleCoverPresentation) = P.n + 1

"""
    numerator_degree(P::DoubleCoverPresentation, m)

Return the degree of the numerator at pole order `m` for `P`.
"""
numerator_degree(P::DoubleCoverPresentation, m::Integer) =
    P.d * m - (P.n + 1 + P.d ÷ 2)

"""
    pole_exponent(P::DoubleCoverPresentation, m)

Return the exponent in the double cover pole `y^(2m-1)`.
"""
pole_exponent(::DoubleCoverPresentation, m::Integer) = 2m - 1

"""
    reduction_scalar(P::DoubleCoverPresentation, m)

Return the exact Griffiths–Dwork scalar for the step from pole order `m + 1`
to `m`.
"""
reduction_scalar(::DoubleCoverPresentation, m::Integer) = QQ(2, 2m - 1)
