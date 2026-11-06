
# (scalar) mul

my_mul!(A::CuModMatrix, B::CuModMatrix, c::Number) = LinearAlgebra.mul!(A, B, c)
my_mul!(A::KaratsubaMatrix, B::KaratsubaMatrix, c::Number) =
    GPUFiniteFieldMatrices.scalar_multiply!(A, B, c)

# matvecmul

my_matvecmul!(z::CuModVector, A::CuModMatrix, b::CuModVector) = LinearAlgebra.mul!(z, A, b)
my_matvecmul!(z::KaratsubaVector, A::KaratsubaMatrix, b::KaratsubaVector) =
    GPUFiniteFieldMatrices.KMatMul!(z, A, b)

# copy

function my_copy!(a, b)
    copy!(a, b)
end

function my_copy!(a::Vector{ZZRingElem}, b::Vector{ZZRingElem})
    for i = 1:length(a)
        Oscar.Nemo.set!(a[i], b[i])
    end
end

function my_copy(a)
    copy(a)
end

function my_copy(a::KaratsubaArray)
    GPUFiniteFieldMatrices.Karatsubacopy(a)
end

function my_copy(a::CuModArray)
    GPUFiniteFieldMatrices.CuModcopy(a)
end

# zero

function my_zero!(a)
    fill!(a, zero(eltype(a)))
    # apparently, the following is not implemented in Julia, it's an Oscar method
    #zero!(a)
end

my_zero!(a::zzModMatrix) = zero!(a)
my_zero!(a::ZZModMatrix) = zero!(a)
my_zero!(a::CuModArray) = GPUFiniteFieldMatrices.zero!(a)
my_zero!(a::KaratsubaArray) = GPUFiniteFieldMatrices.zero!(a)

function my_zero!(a::Vector{ZZRingElem})
    # this is here because of BigInts/ZZRingElem behaving weird in Julia

    for i = 1:length(a)
        Oscar.Nemo.set!(a[i], 0)
    end
end

# sub

function my_sub!(A, B, C)
    Oscar.Nemo.sub!(A, B, C)
end

my_sub!(A::CuModMatrix, B::CuModMatrix, C::CuModMatrix) =
    GPUFiniteFieldMatrices.sub!(A, B, C)
my_sub!(A::KaratsubaArray, B::KaratsubaArray, C::KaratsubaArray) =
    GPUFiniteFieldMatrices.sub!(A, B, C)

function my_sub!(A::ZZModMatrix, B::ZZModMatrix, C::ZZModMatrix)
    @ccall Oscar.Nemo.libflint.fmpz_mod_mat_sub(
        A::Ref{ZZModMatrix},
        B::Ref{ZZModMatrix},
        C::Ref{ZZModMatrix},
        base_ring(B).ninv::Ref{Oscar.Nemo.fmpz_mod_ctx_struct},
    )::Nothing
    return A
end

# add

my_add!(A, B, C) = Oscar.add!(A, B, C)
my_add!(A::CuModArray, B::CuModArray, C::CuModArray) = GPUFiniteFieldMatrices.add!(A, B, C)
my_add!(A::KaratsubaArray, B::KaratsubaArray, C::KaratsubaArray) =
    GPUFiniteFieldMatrices.add!(A, B, C)
