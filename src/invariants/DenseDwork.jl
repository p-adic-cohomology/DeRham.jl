# SPDX-License-Identifier: MIT OR GPL-3.0-or-later
# Integer reconstruction adapted from the workshop's projective_zeta.py.
# All integer reconstruction uses arbitrary precision; see experimental/dense_dwork.md.

function _dense_dwork_input(f)
    K = coefficient_ring(parent(f))
    K isa Oscar.FinField && Oscar.degree(K) == 1 ||
        throw(ArgumentError("dense Dwork requires a prime field"))
    p = BigInt(characteristic(K))
    11 <= p <= 101 || throw(ArgumentError("dense Dwork supports primes 11 through 101"))
    nvars(parent(f)) == 4 && !iszero(f) && is_homogeneous(f) && total_degree(f) == 4 ||
        throw(ArgumentError("dense Dwork requires a homogeneous quartic in four variables"))
    length(f) == 35 || throw(ArgumentError("dense Dwork requires all 35 quartic monomials"))
    return Int(p)
end

function _dense_dwork_matrix(f, p, N)
    executable = joinpath(
        @__DIR__,
        "..",
        "..",
        "deps",
        "usr",
        "bin",
        "dense_dwork" * (Sys.iswindows() ? ".exe" : ""),
    )
    isfile(executable) || error(
        "Dense Dwork backend unavailable; install a C++17 compiler and run Pkg.build(\"DeRham\") once",
    )
    input = IOBuffer()
    println(input, p, " ", N, " 35")
    for i = 1:length(f)
        println(input, join(exponent_vector(f, i), ' '), " ", lift(ZZ, coeff(f, i)))
    end
    seekstart(input)
    output, diagnostic = IOBuffer(), IOBuffer()
    process = run(
        pipeline(
            ignorestatus(`$executable`);
            stdin = input,
            stdout = output,
            stderr = diagnostic,
        ),
    )
    message = String(take!(diagnostic))
    process.exitcode == 3 &&
        throw(ArgumentError("quartic is geometrically singular (Jacobian ranks)"))
    success(process) || error("dense Dwork native computation failed: $message")
    rows = split(strip(String(take!(output))), '\n')
    length(rows) == 43 && rows[1] == "21" || error("invalid dense Dwork matrix response")
    A = [parse.(BigInt, split(rows[i])) for i = 23:43]
    all(row -> length(row) == 21 && all(x -> 0 <= x < big(p)^N, row), A) ||
        error("invalid dense Dwork matrix residues")
    return A
end

function _dense_dwork_residues(A, p, N)
    r = length(A)
    orders = [
        minimum((iszero(A[i][j]) ? N : Int(valuation(ZZ(A[i][j]), p))) for i = 1:r) for
        j = 1:r
    ]
    sort!(orders)
    cp = charpoly(matrix(ZZ, r, r, reduce(vcat, A)))
    residues = Tuple{BigInt,BigInt}[(1, 1)]
    digits = N
    for k = 1:r
        k > 1 && (digits += orders[k-1])
        raw = BigInt(coeff(cp, r-k))
        iszero(mod(raw, big(p)^min(digits, k))) ||
            error("inconsistent scaled Frobenius coefficient at index $k")
        modulus = big(p)^max(0, digits-k)
        push!(residues, (mod(div(raw, big(p)^k), modulus), modulus))
    end
    return residues
end

# Enumerate EVERY prefix allowed by the residues and Newton/Weil bounds.
# A cap requests more precision; it never discards candidates selectively.
function _dense_dwork_lift(residues, p; max_prefixes = 4096)
    r = length(residues)-1
    determinant, determinant_modulus = residues[end]
    signs = [s for s in (1, -1) if iszero(mod(s*big(p)^r-determinant, determinant_modulus))]
    states = [(BigInt[1], BigInt[0])]
    for k = 1:div(r, 2)
        residue, modulus = residues[k+1]
        bound = r*big(p)^k
        next = Tuple{Vector{BigInt},Vector{BigInt}}[]
        for (c, traces) in states
            previous = -sum((c[i+1]*traces[k-i+1] for i = 1:(k-1)); init = big(0))
            lo = cld(previous-bound-k*residue, k*modulus)
            hi = fld(previous+bound-k*residue, k*modulus)
            length(next)+max(0, hi-lo+1) > max_prefixes && return nothing
            for j = lo:hi
                ck = residue+j*modulus
                push!(next, (vcat(c, ck), vcat(traces, previous-k*ck)))
            end
        end
        states = next
    end
    survivors = Vector{BigInt}[]
    for (half, initial_traces) in states, s in signs
        c = vcat(half, [s*big(p)^(r-2i)*half[i+1] for i in reverse(0:div(r, 2))])
        all(iszero(mod(x-a, m)) for (x, (a, m)) in zip(c, residues)) || continue
        traces = copy(initial_traces)
        admissible = true
        for k = length(half):4r
            tr = -sum((c[i+1]*traces[k-i+1] for i = 1:min(r, k-1)); init = big(0))
            k <= r && (tr -= k*c[k+1])
            if abs(tr) > r*big(p)^k
                admissible = false
                break
            end
            push!(traces, tr)
        end
        if admissible
            push!(survivors, c)
            length(survivors) > 1 && return nothing
        end
    end
    isempty(survivors) && error("matrix residues contradict the Weil conditions")
    return only(survivors)
end

"""
    dense_dwork_coefficients(f)

Opt-in experimental dense Dwork–Cartier with Harvey window coefficient extraction.
Return the 22 descending integer coefficients of the degree-21 primitive
Frobenius polynomial of a smooth quartic K3 surface. Requires all 35 monomials
and a prime coefficient field with `11 <= p <= 101`. Singular or unsupported
inputs raise `ArgumentError`. Ambiguity at the native precision limit raises
an error; a candidate is never returned without a unique integer lift.

The optional C++17 backend must have been built once with `Pkg.build("DeRham")`.
Solves need neither a compiler nor Python. See the experimental documentation
for the precision proof, process protocol and supported platforms.
"""
function dense_dwork_coefficients(f::MPolyRingElem)
    p = _dense_dwork_input(f)
    # Julia owns validation. Native Jacobian ranks independently certify smoothness.
    check_smoothness(f) || throw(ArgumentError("quartic is geometrically singular"))
    for N = 3:(p<=47 ? 5 : 4)
        A = _dense_dwork_matrix(f, p, N)
        c = _dense_dwork_lift(_dense_dwork_residues(A, p, N), p)
        isnothing(c) || return ZZ.(reverse(c))
    end
    error("integer lift remains ambiguous at the supported native precision limit")
end

dense_dwork_coefficients(X::ProjectiveHypersurface) =
    dense_dwork_coefficients(polynomial(X))

"""
    dense_dwork_zeta(f)

Return the exact full projective zeta function of a supported dense quartic K3
as an element of `fraction_field(ZZ[T])`, including `(1-T)(1-p*T)(1-p^2*T)`.
Uses [`dense_dwork_coefficients`](@ref); the same restrictions and errors apply.
"""
function dense_dwork_zeta(f::MPolyRingElem)
    cs = dense_dwork_coefficients(f)
    return zeta_function(cs, characteristic(parent(f)), 2)
end

dense_dwork_zeta(X::ProjectiveHypersurface) = dense_dwork_zeta(polynomial(X))
