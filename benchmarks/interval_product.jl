using DeRham
using Oscar

function interval_case(r, modulus)
    ring, = residue_ring(ZZ, modulus)
    matrixspace = matrix_space(ring, r, r)
    a = matrixspace()
    b = matrixspace()
    for j in 1:r, i in 1:r
        a[i, j] = ring(mod(13i + 7j + i * j, 101))
        b[i, j] = ring(mod(11i + 17j + i * j, 101))
    end
    if modulus > ZZ(2)^64
        g = ZZRingElem[ZZ(i * i + 2) for i in 1:r]
    else
        g = UInt[data(ring(i * i + 2)) for i in 1:r]
    end
    return ring, matrixspace, a, b, g
end

function interval_time(method, case, interval_length; repetitions = 3)
    _, matrixspace, a, b, initial = case
    times = Float64[]
    for _ in 1:repetitions
        g = copy(initial)
        temp = matrixspace()
        g_temp = copy(initial)
        GC.gc()
        elapsed = @elapsed if method === :sqrt
            DeRham.sqrt_interval_product!(a, b, 0, interval_length - 1, g, temp, g_temp)
        else
            DeRham.finitediff_prodeval_linear!(
                a,
                b,
                0,
                interval_length - 1,
                g,
                temp,
                g_temp,
            )
        end
        push!(times, elapsed)
    end
    sort!(times)
    return times[cld(length(times), 2)]
end

function fit_power_law(lengths, times)
    x = log.(Float64.(lengths))
    y = log.(times)
    x_mean = sum(x) / length(x)
    y_mean = sum(y) / length(y)
    exponent = sum((x .- x_mean) .* (y .- y_mean)) / sum((x .- x_mean) .^ 2)
    return exponent, y_mean - exponent * x_mean
end

function run_interval_benchmarks()
    lengths = (64, 256, 1024, 4096, 16384)
    moduli = (101^3, ZZ(101)^30, ZZ(2)^80)
    println("r,modulus,L,finitediff_seconds,sqrt_seconds,speedup")
    for r in (2, 4), modulus in moduli
        case = interval_case(r, modulus)
        interval_time(:finitediff, case, 4; repetitions = 1)
        interval_time(:sqrt, case, 4; repetitions = 1)
        finite = Float64[]
        square_root = Float64[]
        for interval_length in lengths
            expected = copy(case[5])
            actual = copy(case[5])
            DeRham.finitediff_prodeval_linear!(
                case[3], case[4], 0, interval_length - 1, expected,
                case[2](), copy(expected),
            )
            DeRham.sqrt_interval_product!(
                case[3], case[4], 0, interval_length - 1, actual,
                case[2](), copy(actual),
            )
            @assert actual == expected
            finite_time = interval_time(:finitediff, case, interval_length)
            sqrt_time = interval_time(:sqrt, case, interval_length)
            push!(finite, finite_time)
            push!(square_root, sqrt_time)
            println(
                "$(r),$(modulus),$(interval_length),$(finite_time),$(sqrt_time),$(finite_time / sqrt_time)",
            )
        end
        crossover = findfirst(i -> square_root[i] < finite[i], eachindex(lengths))
        crossover_length = crossover === nothing ? "none" : string(lengths[crossover])
        finite_exponent, finite_intercept = fit_power_law(lengths[2:end], finite[2:end])
        sqrt_exponent, sqrt_intercept =
            fit_power_law(lengths[2:end], square_root[2:end])
        fitted_crossover = finite_exponent > sqrt_exponent ?
            exp((sqrt_intercept - finite_intercept) / (finite_exponent - sqrt_exponent)) :
            Inf
        println(
            "summary,r=$(r),modulus=$(modulus),crossover_L=$(crossover_length),fitted_crossover_L=$(fitted_crossover),finite_exponent=$(finite_exponent),sqrt_exponent=$(sqrt_exponent)",
        )
    end

    p = 1009
    f = DeRham.fermat_hypersurface(3, 3, p)
    results = Dict{Symbol,Vector{ZZRingElem}}()
    end_to_end = Dict{Symbol,Float64}()
    for method in (:finitediff, :sqrt)
        options = (
            algorithm = :pchunk,
            changef = false,
            interval_product_algorithm = method,
            interval_product_min_length = 1,
        )
        result = DeRham.zeta_coefficients(f; options...)
        results[method] = result
        samples = Float64[]
        for _ in 1:3
            GC.gc()
            push!(samples, @elapsed DeRham.zeta_coefficients(f; options...))
        end
        sort!(samples)
        end_to_end[method] = samples[cld(length(samples), 2)]
    end
    @assert results[:finitediff] == results[:sqrt]
    println(
        "end_to_end,p=$(p),algorithm=pchunk,finitediff_seconds=$(end_to_end[:finitediff]),sqrt_seconds=$(end_to_end[:sqrt]),speedup=$(end_to_end[:finitediff] / end_to_end[:sqrt])",
    )
end

run_interval_benchmarks()
