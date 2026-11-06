
function remove_duplicates_gpu!(costadata_arr)
    i = 1
    while i <= (length(costadata_arr)-1)
        j = i+1
        while j <= length(costadata_arr)
            if all(costadata_arr[i][1][1] .== costadata_arr[j][1][1])
                my_add!(
                    costadata_arr[i][1][2],
                    costadata_arr[i][1][2],
                    costadata_arr[j][1][2],
                )
                deleteat!(costadata_arr, j)
            else
                j = j + 1
            end
        end
        i = i + 1
    end
end

function cuMod(A::zzModMatrix)
    m = modulus(base_ring(parent(A)))
    # if this conversion is illegal,
    #then we're not allowed to make a CuModMatrix
    M = Int(m)

    float_A = float_entries(A)
    CuModMatrix(float_A, M, elem_type = Float64)
end

function KaratsubaMat(A::zzModMatrix, A_temp = nothing)
    m = modulus(base_ring(parent(A)))
    temp = collect(factor(m))[1]
    d = temp[2]
    p = temp[1]

    N1 = Int(p)^Int(round(d/2))
    N2 = Int(p)^Int(d - round(d/2))

    GPUFiniteFieldMatrices.KaratsubaMatrix(Float64, cuMod(A), N1, N2, N1*N2)
end

function Karatsuba_copyto!(A_gpu::KaratsubaMatrix, A::Matrix{Float64}, A_temp = nothing)

    temp = CuModMatrix(A, A_gpu.N1*A_gpu.N2, elem_type = Float64)

    GPUFiniteFieldMatrices.divide_elements!(A_gpu.data2, temp, A_gpu.N1)
    GPUFiniteFieldMatrices.mod_elements!(A_gpu.data2, A_gpu.N2)

    LinearAlgebra.mul!(A_gpu.data1, A_gpu.data2, Int(A_gpu.N1))
    GPUFiniteFieldMatrices.sub!(A_gpu.data1, temp, A_gpu.data1; mod_N = A_gpu.N1)
    GPUFiniteFieldMatrices.mod_elements!(A_gpu.data1, A_gpu.N1)

    A_gpu
end

#TODO: incorporate this idiomatically using adapt.jl
function adapt_to_gpu(Ruvs)
    gpu_Ruvs = Dict{Vector{Int64},Vector{CuModMatrix{Float64}}}()

    for (v, linpoly) in Ruvs

        gpu_Ruvs[v] = Vector(undef, length(linpoly))
        for i = 1:length(linpoly)
            gpu_Ruvs[v][i] = cuMod(linpoly[i])
        end
    end

    gpu_Ruvs
end


"""
Return a Matrix{Float64} which has entries that are the same
as A
"""
function float_entries(A::zzModMatrix)
    res = zeros(Float64, size(A)...)

    for i = 1:size(A, 1)
        for j = 1:size(A, 2)
            dat = data(A[i, j])
            # lifted = lift(ZZ,A[i,j])
            res[i, j] = convert(Float64, dat)
        end
    end

    res
end

"""
Updates dest with the entries of src as Float64 values
"""
function float_entries!(dest::Matrix{Float64}, src::zzModMatrix)
    res = dest
    A = src

    for i = 1:size(A, 1)
        for j = 1:size(A, 2)
            lifted = lift(ZZ, A[i, j])
            res[i, j] = convert.(Float64, convert.(Int64, lifted))
        end
    end

    res
end

"""
Updates dest with the entries of src as Float64 values
"""
function float_entries!(dest::Matrix{Float64}, src::ZZMatrix)
    res = dest
    A = src

    for i = 1:size(A, 1)
        for j = 1:size(A, 2)
            res[i, j] = convert(Float64, convert(Int64, A[i, j]))
        end
    end

    res
end
