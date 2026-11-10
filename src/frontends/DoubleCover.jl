"""
    DoubleCover(F)

The zeta frontend for the double cover of projective space branched along `F`.
"""
struct DoubleCover{F} <: ZetaFrontend
    F::F
    P::DoubleCoverPresentation
end

function DoubleCover(F)
    n = nvars(parent(F)) - 1
    d = total_degree(F)
    return DoubleCover{typeof(F)}(F, DoubleCoverPresentation(n, d))
end

DoubleCover{F}(f::F) where {F} = DoubleCover(f)

polynomial(fe::DoubleCover) = fe.F
presentation(fe::DoubleCover) = fe.P
weight(fe::DoubleCover) = fe.P.n

function check_input(fe::DoubleCover)
    F = polynomial(fe)
    iseven(total_degree(F)) || throw(ArgumentError("the branch degree must be even"))
    is_homogeneous(F) || throw(ArgumentError("the branch polynomial must be homogeneous"))
    characteristic(parent(F)) == 2 &&
        throw(ArgumentError("double covers are not supported in characteristic 2"))
    return nothing
end

function reduction_prescale(fe::DoubleCover, m::Integer)
    m >= 1 || throw(ArgumentError("the pole order must be positive"))
    scale = QQ(1)
    P = presentation(fe)
    for j in 1:(m-1)
        scale *= inv(reduction_scalar(P, j))
    end
    return scale
end

frobenius_series(fe::DoubleCover, N, m) =
    computeD_half_integral(N, pole_exponent(presentation(fe), m))

function frobenius_term_pole_order(fe::DoubleCover, m, j)
    p = Int(characteristic(parent(polynomial(fe))))
    return (p * (2m - 1 + 2j) + 1) ÷ 2
end
