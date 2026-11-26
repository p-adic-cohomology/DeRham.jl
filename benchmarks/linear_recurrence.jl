using DeRham
using Oscar

function recurrence_case(r, modulus)
    ring, = residue_ring(ZZ, modulus)
    matrixspace = matrix_space(ring, r, r)
    a = matrixspace()
    b = matrixspace()
    for j = 1:r, i = 1:r
        a[i, j] = ring(mod(13i + 7j + i * j, 101))
        b[i, j] = ring(mod(11i + 17j + i * j, 101))
    end
    if modulus > ZZ(2)^64
        g = ZZRingElem[ZZ(i * i + 2) for i = 1:r]
    else
        g = UInt[data(ring(i * i + 2)) for i = 1:r]
    end
    return ring, matrixspace, a, b, g
end

function linear_recurrence_time(method, case, run_length; repetitions = 3)
    _, matrixspace, a, b, initial = case
    times = Float64[]
    for _ = 1:repetitions
        g = deepcopy(initial)
        temp = matrixspace()
        g_temp = deepcopy(initial)
        GC.gc()
        elapsed = @elapsed if method === :sqrt
            DeRham.sqrt_linear_recurrence!(a, b, 0, run_length - 1, g, temp, g_temp)
        else
            DeRham.finitediff_prodeval_linear!(a, b, 0, run_length - 1, g, temp, g_temp)
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

function run_linear_recurrence_benchmarks()
    lengths = (64, 256, 1024, 4096, 16384, 65536, 262144, 1048576)
    moduli = (101^3, ZZ(101)^30, ZZ(2)^80)
    println("r,modulus,L,finitediff_seconds,sqrt_seconds,speedup")
    for r in (2, 4, 10, 28)
        benchmark_moduli = r <= 4 ? moduli : (moduli[1],)
        repetitions = r <= 4 ? 3 : 1
        for modulus in benchmark_moduli
            case = recurrence_case(r, modulus)
            linear_recurrence_time(:finitediff, case, 4; repetitions = 1)
            linear_recurrence_time(:sqrt, case, 4; repetitions = 1)
            finite = Float64[]
            square_root = Float64[]
            measured_lengths = Int[]
            for run_length in lengths
                expected = deepcopy(case[5])
                actual = deepcopy(case[5])
                DeRham.finitediff_prodeval_linear!(
                    case[3],
                    case[4],
                    0,
                    run_length - 1,
                    expected,
                    case[2](),
                    deepcopy(expected),
                )
                DeRham.sqrt_linear_recurrence!(
                    case[3],
                    case[4],
                    0,
                    run_length - 1,
                    actual,
                    case[2](),
                    deepcopy(actual),
                )
                @assert actual == expected
                finite_time = linear_recurrence_time(
                    :finitediff,
                    case,
                    run_length;
                    repetitions = repetitions,
                )
                sqrt_time = linear_recurrence_time(
                    :sqrt,
                    case,
                    run_length;
                    repetitions = repetitions,
                )
                push!(finite, finite_time)
                push!(square_root, sqrt_time)
                push!(measured_lengths, run_length)
                println(
                    "$(r),$(modulus),$(run_length),$(finite_time),$(sqrt_time),$(finite_time / sqrt_time)",
                )
                max(finite_time, sqrt_time) > 180 && break
            end
            crossover =
                findfirst(i -> square_root[i] < finite[i], eachindex(measured_lengths))
            crossover_length =
                crossover === nothing ? "none" : string(measured_lengths[crossover])
            if length(measured_lengths) >= 3
                finite_exponent, finite_intercept =
                    fit_power_law(measured_lengths[2:end], finite[2:end])
                sqrt_exponent, sqrt_intercept =
                    fit_power_law(measured_lengths[2:end], square_root[2:end])
                fitted_crossover =
                    finite_exponent > sqrt_exponent ?
                    exp(
                        (sqrt_intercept - finite_intercept) /
                        (finite_exponent - sqrt_exponent),
                    ) : Inf
            else
                finite_exponent = NaN
                sqrt_exponent = NaN
                fitted_crossover = NaN
            end
            println(
                "summary,r=$(r),modulus=$(modulus),tested_lengths=$(join(measured_lengths, ';')),crossover_L=$(crossover_length),fitted_crossover_L=$(fitted_crossover),finite_exponent=$(finite_exponent),sqrt_exponent=$(sqrt_exponent)",
            )
        end
    end

    p = 1009
    f = DeRham.fermat_hypersurface(3, 3, p)
    results = Dict{Symbol,Vector{ZZRingElem}}()
    end_to_end = Dict{Symbol,Float64}()
    for method in (:finitediff, :sqrt)
        options = (
            algorithm = :pchunk,
            changef = false,
            linear_recurrence_algorithm = method,
            sqrt_linear_recurrence_min_length = 1,
        )
        result = DeRham.zeta_coefficients(f; options...)
        results[method] = result
        samples = Float64[]
        for _ = 1:3
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

run_linear_recurrence_benchmarks()
