"""
    ProjectiveHypersurface(f)

The zeta frontend for the hypersurface defined by `f`.
"""
struct ProjectiveHypersurface{F} <: ZetaFrontend
    f::F
    P::HypersurfaceQuotient
end

function ProjectiveHypersurface(f)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    return ProjectiveHypersurface{typeof(f)}(f, HypersurfaceQuotient(n, d))
end

ProjectiveHypersurface{F}(f::F) where {F} = ProjectiveHypersurface(f)

polynomial(fe::ProjectiveHypersurface) = fe.f
cohomology_quotient(fe::ProjectiveHypersurface) = fe.P
weight(fe::ProjectiveHypersurface) = fe.P.n - 1

check_input(::ProjectiveHypersurface) = nothing

_reduction_allowed_policies(::ProjectiveHypersurface) =
    (:depthfirst, :pchunk, :varbyvar, :akr, :akr_naive, :akr_levelbylevel)

function varbyvar_reduction_setup(fe::ProjectiveHypersurface, S)
    f = polynomial(fe)
    P = cohomology_quotient(fe)
    n = P.n
    required_S = [n]
    if S == required_S
        return :varbyvar, S, nothing, required_S
    elseif n > 4
        for i = 0:n
            if is_Ssmooth(f, P, [n-i])
                return :varbyvar, required_S, i, required_S
            elseif i == n
                return :depthfirst, S, nothing, required_S
            end
        end
    end
    return :depthfirst, S, nothing, required_S
end

reduction_prescale(fe::ProjectiveHypersurface, m::Integer) = reduction_prescale(fe.P, m)

frobenius_series(::ProjectiveHypersurface, N, m) = computeD(N, m)

frobenius_term_pole_order(fe::ProjectiveHypersurface, m, j) =
    Int(characteristic(parent(polynomial(fe)))) * (m + j)
