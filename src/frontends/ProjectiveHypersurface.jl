"""
    ProjectiveHypersurface(f)

The zeta frontend for the hypersurface defined by `f`.
"""
struct ProjectiveHypersurface{F} <: ZetaFrontend
    f::F
    P::HypersurfacePresentation
end

function ProjectiveHypersurface(f)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    return ProjectiveHypersurface{typeof(f)}(f, HypersurfacePresentation(n, d))
end

ProjectiveHypersurface{F}(f::F) where {F} = ProjectiveHypersurface(f)

polynomial(fe::ProjectiveHypersurface) = fe.f
presentation(fe::ProjectiveHypersurface) = fe.P
weight(fe::ProjectiveHypersurface) = fe.P.n - 1

check_input(::ProjectiveHypersurface) = nothing

reduction_prescale(::ProjectiveHypersurface, m::Integer) = factorial(ZZ(m - 1))

frobenius_series(::ProjectiveHypersurface, N, m) = computeD(N, m)

frobenius_term_pole_order(fe::ProjectiveHypersurface, m, j) =
    characteristic(parent(polynomial(fe))) * (m + j)
