
"""
Gives a special PEP for when S = [0] and d=3

This should be all the Ruv ever needed in this situation.
"""
function cubic_S_zero_Vs(n)
    Vs = Vector{Vector{Int}}()

    for i = 1:(n+1)
        V = zeros(Int, n+1)
        V[i] = 3
        push!(Vs, V)
    end

    for i = 1:n
        V1 = zeros(Int, n+1)
        V1[i] = 1
        V1[i+1] = 2
        push!(Vs, V1)

        V2 = zeros(Int, n+1)
        V2[i] = 2
        V2[i+1] = 1
        push!(Vs, V2)
    end

    Vs
end

"""
This function exists so that we can have the option to time it
in verbose mode.
"""
function select_Ruv_PEP(n, d, S, params, compute, lazy, oscar_matspace, cache)

    m = Integer(modulus(base_ring(oscar_matspace)))

    compute_float = V -> float_entries.(compute(V))

    function compute_gpu(V; copyto = nothing)

        if (3 < params.verbose)
            println("Moving Ruv for $V to the gpu")
        end

        if (ZZ(2)^25 < m < ZZ(2)^106)
            if copyto != nothing
                Karatsuba_copyto!.(copyto, float_entries.(compute(V)))
            else
                KaratsubaMat.(compute(V))
            end
        else
            if copyto != nothing
                copyto!.(copyto, float_entries.(compute(V)))
                copyto
            else
                cuMod.(compute(V))
            end
        end
    end


    EXTRA_MEMORY = false

    if 5 == n &&
       d == 3 &&
       S == [n] &&
       params.use_gpu &&
       (ZZ(2)^25 < m < ZZ(2)^106) &&
       !EXTRA_MEMORY
        println("hey howdy hey")

        Ruv = LRULazyPEP{KaratsubaMatrix{Float64}}(compute_gpu, d, usethreads = false)
    elseif 5 == n && d == 3 && S == [n] && params.use_gpu && !EXTRA_MEMORY
        println("howdy")

        Ruv = LRULazyPEP{CuModMatrix{Float64}}(compute_gpu, d, usethreads = false)

    elseif 3 < n && d == 3 && S == [n] && params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106) #4 < n
        println("hey")

        eager_Vs = cubic_S_zero_Vs(n)

        cpu_Ruv =
            LazyPEP{Matrix{Float64}}(compute_float, eagerVs = eager_Vs, usethreads = false)
        m = M = Int(modulus(base_ring(oscar_matspace)))
        temp = collect(factor(m))[1]
        d = temp[2]
        p = temp[1]

        N1 = Int(p)^Int(round(d/2))
        N2 = Int(p)^Int(d - round(d/2))

        if (3 < params.verbose)
            create_gpu =
                matrices -> begin
                    println("Moving an Ruv to the GPU (first time creation!)")
                    CUDA.@time KaratsubaMatrix.(
                        Float64,
                        CuModMatrix.(matrices, N1*N2, elem_type = Float64),
                        N1,
                        N2,
                        N1*N2,
                    )
                end

            convert_gpu =
                (dest, matrices) -> begin
                    println("Moving an Ruv to the GPU (already allocated)")

                    CUDA.@time begin
                        for i = 1:length(dest)
                            Karatsuba_copyto!(dest[i], matrices[i])
                        end
                    end
                end
        else
            # should it really be N1??? or actually N1 * N2 TODO
            create_gpu =
                matrices -> KaratsubaMatrix.(
                    Float64,
                    CuModMatrix.(matrices, N1*N2, elem_type = Float64),
                    N1,
                    N2,
                    N1*N2,
                )
            convert_gpu =
                (dest, matrices) -> begin
                    for i = 1:length(dest)
                        Karatsuba_copyto!(dest[i], matrices[i])
                    end
                end
        end
        s = size(oscar_matspace(), 1)

        # memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
        # # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
        # memory = s^2 * 8 * (n+2)

        # maxsize = div(memory_cap,memory)

        #println(maxsize)
        #testing
        if n == 5 # (cubic) fourfold
            maxsize = 3
        elseif n == 4 #(cubic) threefold
            maxsize = 3#20
        else
            memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
            # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
            memory = s^2 * 8 * (n+2)

            maxsize = div(memory_cap, memory)
        end

        Ruv = CachePEP{Matrix{Float64},KaratsubaMatrix{Float64}}(
            cpu_Ruv,
            create_gpu,
            convert_gpu,
            maxsize,
        )

    elseif 3 < n && d == 3 && S == [n] && params.use_gpu #4 < n
        println("hello")

        eager_Vs = cubic_S_zero_Vs(n)

        cpu_Ruv =
            LazyPEP{Matrix{Float64}}(compute_float, eagerVs = eager_Vs, usethreads = false)
        m = M = Int(modulus(base_ring(oscar_matspace)))

        if (3 < params.verbose)
            create_gpu =
                matrices -> begin
                    println("Moving an Ruv to the GPU (first time creation!)")
                    CUDA.@time CuModMatrix.(matrices, M, elem_type = Float64)
                end

            convert_gpu =
                (dest, matrices) -> begin
                    println("Moving an Ruv to the GPU (already allocated)")

                    CUDA.@time begin
                        for i = 1:length(dest)
                            copyto!(dest[i], matrices[i])
                        end
                    end
                end
        else
            create_gpu = matrices -> CuModMatrix.(matrices, M, elem_type = Float64)
            convert_gpu = (dest, matrices) -> begin
                for i = 1:length(dest)
                    copyto!(dest[i], matrices[i])
                end
            end
        end
        s = size(oscar_matspace(), 1)

        if n == 5 # (cubic) fourfold
            maxsize = 3 # really should probably make this an LRU cache for varbyvar
        elseif n == 4 #(cubic) threefold
            maxsize = 20
        else
            memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
            # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
            memory = s^2 * 8 * (n+2)

            maxsize = div(memory_cap, memory)
        end
        #maxsize = 13

        Ruv = CachePEP{Matrix{Float64},CuModMatrix{Float64}}(
            cpu_Ruv,
            create_gpu,
            convert_gpu,
            maxsize,
        )

        #(1 < params.verbose) && println("Initial cache info: $(cache_info(Ruv.Ucomponent))")
    elseif params.use_gpu && lazy && (ZZ(2)^25 < m < ZZ(2)^106)
        compute_gpu_karatsuba = V -> @. KaratsubaMat(float_entries(compute(V)))
        Ruv = LazyPEP{KaratsubaMatrix{Float64}}(compute_gpu_karatsuba)

    elseif params.use_gpu && lazy
        Ruv = LazyPEP{CuModMatrix{Float64}}(compute_gpu)

    elseif params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106)
        compute_gpu_karatsuba = V -> @. KaratsubaMat(float_entries(compute(V)))
        Ruv = EagerPEP{KaratsubaMatrix{Float64}}(
            cache[d],
            compute_gpu_karatsuba,
            usethreads = false,
        )

    elseif params.use_gpu
        Ruv = EagerPEP{CuModMatrix{Float64}}(
            cache[d],
            compute_gpu,
            usethreads = params.use_threads,
        )

    elseif lazy
        Ruv = LazyPEP{typeof(oscar_matspace())}(compute)

    else
        Ruv = EagerPEP{typeof(oscar_matspace())}(
            cache[d],
            compute,
            usethreads = params.use_threads,
        )
    end

    Ruv
end

function pregen_select_Ruv_PEP(n, d, S, params, oscar_matspace)

    m = Integer(modulus(base_ring(oscar_matspace)))

    if 3 < n && d == 3 && S == [n] && params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106) #4 < n

        cpu_Ruv = PregenLazyPEP{Matrix{Float64}}(nothing)
        m = M = Int(modulus(base_ring(oscar_matspace)))
        temp = collect(factor(m))[1]
        d = temp[2]
        p = temp[1]

        N1 = Int(p)^Int(round(d/2))
        N2 = Int(p)^Int(d - round(d/2))

        create_gpu =
            matrices -> KaratsubaMatrix.(
                Float64,
                CuModMatrix.(matrices, N1, elem_type = Float64),
                N1,
                N2,
                N1*N2,
            )
        convert_gpu = (dest, matrices) -> begin
            for i = 1:length(dest)
                copyto!(dest[i], matrices[i])
            end
        end
        s = size(oscar_matspace(), 1)

        memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
        # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
        memory = s^2 * 8 * (n+2)

        maxsize = div(memory_cap, memory)

        #println(maxsize)
        #testing
        #maxsize = 13

        Ruv = CachePEP{Matrix{Float64},KaratsubaMatrix{Float64}}(
            cpu_Ruv,
            create_gpu,
            convert_gpu,
            maxsize,
        )

    elseif 3 < n && d == 3 && S == [n] && params.use_gpu #4 < n

        cpu_Ruv = PregenLazyPEP{Matrix{Float64}}(nothing)
        m = M = Int(modulus(base_ring(oscar_matspace)))

        create_gpu = matrices -> CuModMatrix.(matrices, M, elem_type = Float64)
        convert_gpu = (dest, matrices) -> begin
            for i = 1:length(dest)
                copyto!(dest[i], matrices[i])
            end
        end
        s = size(oscar_matspace(), 1)

        memory_cap = totalmem(CUDA.device())#4_500_000_000 # about 6 gigabytes of gpu memory
        # 8 bytes per float64, n+2 matrices, s^2 entries per matrix
        memory = s^2 * 8 * (n+2)

        maxsize = div(memory_cap, memory)

        #println(maxsize)
        #testing
        #maxsize = 13

        Ruv = CachePEP{Matrix{Float64},CuModMatrix{Float64}}(
            cpu_Ruv,
            create_gpu,
            convert_gpu,
            maxsize,
        )

        #(1 < params.verbose) && println("Initial cache info: $(cache_info(Ruv.Ucomponent))")

    elseif params.use_gpu && (ZZ(2)^25 < m < ZZ(2)^106)
        Ruv = PregenLazyPEP{KaratsubaMatrix{Float64}}(nothing)
    elseif params.use_gpu
        Ruv = PregenLazyPEP{CuModMatrix{Float64}}(nothing)
    else
        Ruv = PregenLazyPEP{typeof(oscar_matspace())}(nothing)
    end

    Ruv
end
