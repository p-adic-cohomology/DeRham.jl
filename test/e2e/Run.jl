struct CaseOutcome
    name::String
    passed::Bool
    skipped::Bool
    skip_reason::Union{Nothing,String}
    expected_runtime::Float64
    actual_runtime::Float64
end

function runtime_drift(expected::Real, actual::Real)
    expected == 0 && return actual > 0
    return actual > 2 * expected || actual < expected / 2
end

runtime_drift(outcome::CaseOutcome) =
    !outcome.skipped && runtime_drift(outcome.expected_runtime, outcome.actual_runtime)

function print_runtime_summary(outcomes::AbstractVector{CaseOutcome}, io::IO = stdout)
    println(io, "[e2e] runtime summary")
    @printf(io, "%-52s %12s %12s %s\n", "record", "declared (s)", "measured (s)", "flag",)
    for outcome in outcomes
        flag = if outcome.skipped
            "SKIP"
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

function run_selected(
    catalogue::ExampleCatalogue,
    selected::NamedTuple,
    expected::Dict{String,Any},
)
    records = selected.records
    broad = selected.broad

    estimated_total = isempty(records) ? 0.0 : sum(r.runtime for r in records)
    println(
        "[e2e] selected $(length(records)) case(s), estimated total runtime $(round(estimated_total, digits = 3))s",
    )

    outcomes = CaseOutcome[]
    for record in records
        caps = Common.required_capabilities(record)
        unavailable = [c for c in caps if !Common.capability_available(c)]

        if !isempty(unavailable)
            reason = "requires unavailable capabilities: $(sort(collect(unavailable)))"
            if broad
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
        executor.execute(F, record.parameters)
        actual_runtime = @elapsed raw = executor.execute(F, record.parameters)
        actual = executor.normalize(raw)
        passed = executor.compare(expected[record.name], actual)
        status = passed ? "PASS" : "FAIL"
        println(
            "[e2e] $(status) $(record.name)  expected_runtime=$(record.runtime)s actual_runtime=$(round(actual_runtime, digits = 3))s",
        )
        @test passed
        push!(
            outcomes,
            CaseOutcome(
                record.name,
                passed,
                false,
                nothing,
                record.runtime,
                actual_runtime,
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
    return run_selected(catalogue, selected, expected)
end
