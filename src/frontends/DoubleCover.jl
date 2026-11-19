"""
    DoubleCover(F)

The zeta frontend for the double cover of projective space branched along `F`.
"""
struct DoubleCover{F} <: ZetaFrontend
    F::F
    P::DoubleCoverQuotient
end

function DoubleCover(F)
    n = nvars(parent(F)) - 1
    d = total_degree(F)
    return DoubleCover{typeof(F)}(F, DoubleCoverQuotient(n, d))
end

DoubleCover{F}(f::F) where {F} = DoubleCover(f)

polynomial(fe::DoubleCover) = fe.F
cohomology_quotient(fe::DoubleCover) = fe.P
weight(fe::DoubleCover) = fe.P.n

function check_input(fe::DoubleCover)
    F = polynomial(fe)
    iseven(total_degree(F)) || throw(ArgumentError("the branch degree must be even"))
    is_homogeneous(F) || throw(ArgumentError("the branch polynomial must be homogeneous"))
    p = characteristic(parent(F))
    p == 0 && throw(
        ArgumentError(
            "double covers must be defined over a finite field of odd characteristic",
        ),
    )
    p == 2 && throw(ArgumentError("double covers are not supported in characteristic 2"))
    return nothing
end

_reduction_allowed_policies(::DoubleCover) = (:depthfirst, :varbyvar)

reduction_prescale(fe::DoubleCover, m::Integer) = reduction_prescale(fe.P, m)

frobenius_series(fe::DoubleCover, N, m) =
    computeD_half_integral(N, pole_exponent(cohomology_quotient(fe), m))

function frobenius_term_pole_order(fe::DoubleCover, m, j)
    p = Int(characteristic(parent(polynomial(fe))))
    return (p * (2m - 1 + 2j) + 1) ÷ 2
end
