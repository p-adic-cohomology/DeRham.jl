struct CaseOutcome
    name::String
    passed::Bool
    skipped::Bool
    skip_reason::Union{Nothing,String}
    error_message::Union{Nothing,String}
    expected_runtime::Float64
    actual_runtime::Float64
end

CaseOutcome(
    name::AbstractString,
    passed::Bool,
    skipped::Bool,
    skip_reason::Union{Nothing,AbstractString},
    expected_runtime::Real,
    actual_runtime::Real,
) = CaseOutcome(
    String(name),
    passed,
    skipped,
    skip_reason === nothing ? nothing : String(skip_reason),
    nothing,
    Float64(expected_runtime),
    Float64(actual_runtime),
)

function runtime_drift(expected::Real, actual::Real)
    expected == 0 && return actual > 0
    return actual > 2 * expected || actual < expected / 2
end

runtime_drift(outcome::CaseOutcome) =
    !outcome.skipped &&
    outcome.error_message === nothing &&
    runtime_drift(outcome.expected_runtime, outcome.actual_runtime)

function print_runtime_summary(outcomes::AbstractVector{CaseOutcome}, io::IO = stdout)
    println(io, "[e2e] runtime summary")
    @printf(io, "%-52s %12s %12s %s\n", "record", "declared (s)", "measured (s)", "flag",)
    for outcome in outcomes
        flag = if outcome.skipped
            "SKIP"
        elseif outcome.error_message !== nothing
            "ERROR"
        elseif !outcome.passed
            "FAIL"
        elseif runtime_drift(outcome)
            "DRIFT"
        else
            ""
        end
        @printf(
            io,
            "%-52s %12.3f %12.3f %s\n",
            outcome.name,
            outcome.expected_runtime,
            outcome.actual_runtime,
            flag,
        )
    end
    return nothing
end

function _capture_case_step(stage::AbstractString, f::Function)
    try
        return (value = f(), error_message = nothing)
    catch exception
        exception isa InterruptException && rethrow()
        return (
            value = nothing,
            error_message = "$(stage): $(sprint(showerror, exception))",
        )
    end
end

function _timed_execute(executor, input, parameters::AbstractDict, stage::AbstractString)
    started_at = time_ns()
    try
        value = executor.execute(input, parameters)
        return (
            value = value,
            error_message = nothing,
            runtime = (time_ns() - started_at) / 1.0e9,
        )
    catch exception
        exception isa InterruptException && rethrow()
        runtime = (time_ns() - started_at) / 1.0e9
        return (
            value = nothing,
            error_message = "$(stage): $(sprint(showerror, exception))",
            runtime = runtime,
        )
    end
end

function _case_error(error_message, actual_runtime)
    return (
        passed = false,
        errored = true,
        error_message = error_message,
        actual_runtime = actual_runtime,
    )
end

function _comparison_result(executor, expected, actual, stage::AbstractString)
    comparison = _capture_case_step(stage, () -> executor.compare(expected, actual))
    comparison.error_message === nothing || return comparison
    comparison.value isa Bool || return (
        value = nothing,
        error_message = "$(stage): comparison returned $(typeof(comparison.value)); expected Bool",
    )
    return comparison
end

function evaluate_case(executor, input, parameters::AbstractDict, expected_value)
    execution = _timed_execute(executor, input, parameters, "execution")
    execution.error_message === nothing ||
        return _case_error(execution.error_message, execution.runtime)

    normalization =
        _capture_case_step("normalization", () -> executor.normalize(execution.value))
    normalization.error_message === nothing ||
        return _case_error(normalization.error_message, execution.runtime)

    comparison =
        _comparison_result(executor, expected_value, normalization.value, "comparison")
    comparison.error_message === nothing ||
        return _case_error(comparison.error_message, execution.runtime)

    return (
        passed = comparison.value,
        errored = false,
        error_message = nothing,
        result = normalization.value,
        actual_runtime = execution.runtime,
    )
end

function run_selected(
    catalogue::ExampleCatalogue,
    selected::NamedTuple,
    expected::Dict{String,Any};
    priming_records::Vector{ExecutionRecord} = ExecutionRecord[],
    capability_available::Function = Common.capability_available,
)
    primed_names = Set(record.name for record in priming_records)
    selected_names = Set(record.name for record in selected.records)
    records = vcat(
        priming_records,
        filter(record -> !(record.name in primed_names), selected.records),
    )

    estimated_total = isempty(records) ? 0.0 : sum(r.runtime for r in records)
    isempty(priming_records) ||
        println("[e2e] priming set: $(length(priming_records)) record(s)")
    println(
        "[e2e] selected $(length(selected.records)) case(s), running $(length(records)) record(s), estimated total runtime $(round(estimated_total, digits = 3))s",
    )

    outcomes = CaseOutcome[]
    for record in records
        caps = Common.required_capabilities(record)
        unavailable = [c for c in caps if !capability_available(c)]

        if !isempty(unavailable)
            reason = "requires unavailable capabilities: $(sort(collect(unavailable)))"
            if selected.broad ||
               (record.name in primed_names && !(record.name in selected_names))
                println("[e2e] SKIP $(record.name): $(reason)")
                push!(
                    outcomes,
                    CaseOutcome(record.name, false, true, reason, record.runtime, 0.0),
                )
                continue
            else
                error("execution record $(repr(record.name)) $(reason)")
            end
        end

        example = catalogue.examples[record.example]
        executor = EXECUTOR_REGISTRY[record.type]
        F = Common.build_input(example, record)
        result = evaluate_case(executor, F, record.parameters, expected[record.name])
        status = result.errored ? "ERROR" : result.passed ? "PASS" : "FAIL"
        line =
            "[e2e] $(status) $(record.name)  " *
            "expected_runtime=$(record.runtime)s " *
            "actual_runtime=$(round(result.actual_runtime, digits = 3))s"
        if result.errored
            line *= "  $(result.error_message)"
        end
        println(line)
        if result.errored
            @test "E2E record $(record.name): $(result.error_message)" == ""
        else
            @test result.passed
        end
        push!(
            outcomes,
            CaseOutcome(
                record.name,
                result.passed,
                false,
                nothing,
                result.error_message,
                record.runtime,
                result.actual_runtime,
            ),
        )
    end
    print_runtime_summary(outcomes)
    return outcomes
end

"""
Runs the selected E2E records and prints the runtime summary table.
"""
function run_e2e(e2e_root::AbstractString, args::AbstractVector{<:AbstractString})
    config = load_config(joinpath(e2e_root, "config.toml"))
    catalogue = load_example_catalogue(e2e_root)
    records = load_execution_records(e2e_root)
    expected = validate_catalogue(catalogue, records, config)
    selection = parse_selection(args)
    selected = select_records(records, config, selection)
    priming_records = if should_run_priming(selection)
        records_by_name = Dict(record.name => record for record in records)
        [records_by_name[name] for name in config.priming_records]
    else
        ExecutionRecord[]
    end
    return run_selected(catalogue, selected, expected; priming_records = priming_records)
end
