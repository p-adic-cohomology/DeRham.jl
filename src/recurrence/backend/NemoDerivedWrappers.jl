# SPDX-License-Identifier: BSD-2-Clause

# Derived from Nemo.jl (https://github.com/Nemocas/Nemo.jl, src/flint/nmod_mat.jl and src/flint/fmpz_mod_mat.jl), Copyright (c) 2014-2024: William Hart, Tommy Hofmann, Claus Fieker, Fredrik Johansson, and Nemo developers.

# Copyright (c) 2025 the DeRham.jl contributors. See LICENSES/BSD-2-Clause.txt.

function my_addmul!(A::zzModMatrix, B::zzModMatrix, C::UInt, D::zzModMatrix)
    @ccall Oscar.Nemo.libflint.nmod_mat_scalar_addmul_ui(
        A::Ref{zzModMatrix},
        B::Ref{zzModMatrix},
        D::Ref{zzModMatrix},
        C::UInt,
    )::Nothing
    return A
end

function my_mul!(A::zzModMatrix, B::zzModMatrix, c)
    n = characteristic(base_ring(parent(A)))
    ui(i) = 0 ≤ i ? UInt(i % n) : UInt((i % n) + n)
    c = ui(c)
    @ccall Oscar.Nemo.libflint.nmod_mat_scalar_mul(
        A::Ref{zzModMatrix},
        B::Ref{zzModMatrix},
        c::UInt,
    )::Nothing
    return A
end

function my_mul!(A::ZZModMatrix, B::ZZModMatrix, c::Int)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_scalar_mul_si(
        A::Ref{ZZModMatrix},
        B::Ref{ZZModMatrix},
        c::Int,
        base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
    )::Nothing
    return A
end

function my_mul!(A::ZZModMatrix, B::ZZModMatrix, c::ZZRingElem)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_scalar_mul_fmpz(
        A::Ref{ZZModMatrix},
        B::Ref{ZZModMatrix},
        C::Ref{ZZRingElem},
        base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
    )::Nothing
    return A
end

function my_matvecmul!(z::Vector{UInt}, A::zzModMatrix, b::Vector{UInt})
    @ccall Oscar.Nemo.libflint.nmod_mat_mul_nmod_vec(
        z::Ptr{UInt},
        A::Ref{zzModMatrix},
        b::Ptr{UInt},
        length(b)::Int,
    )::Nothing
    return z
end

function my_matvecmul!(z::Vector{ZZRingElem}, A::ZZModMatrix, b::Vector{ZZRingElem})
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_mul_fmpz_vec_ptr(
        z::Ptr{Ref{ZZRingElem}},
        A::Ref{ZZModMatrix},
        b::Ptr{Ref{ZZRingElem}},
        length(b)::Int,
        base_ring(A).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
    )::Nothing
    return z
end
