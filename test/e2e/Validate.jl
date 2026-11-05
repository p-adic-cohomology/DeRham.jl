# Validate.jl — validates the entire combined example + execution-record
# catalogue before any selected case runs. Cheap: no Oscar objects are built
# and no zeta computation happens here (see executors/ZetaCoefficients.jl).
# Errors are aggregated across the whole catalogue and reported together,
# each tagged with its source file and record/example id, rather than
# failing on the first problem found.

struct CatalogueValidationError <: Exception
    errors::Vector{String}
end

function Base.showerror(io::IO, e::CatalogueValidationError)
    println(io, "E2E catalogue validation failed with $(length(e.errors)) error(s):")
    for msg in e.errors
        println(io, "  - $(msg)")
    end
end

const PRIMING_CODE_PATH_DEFAULTS = (
    "algorithm" => "default",
    "fastevaluation" => true,
    "changef" => true,
    "always_use_bigints" => false,
    "use_gpu" => false,
    "use_threads" => false,
    "termorder" => "invlex",
    "k" => nothing,
)

function _canonical_priming_option(key::String, value)
    if key in ("algorithm", "termorder") && value isa AbstractString
        return String(value)
    elseif key == "k" && value isa Integer
        # `k` changes the iteration bound; its integer type determines specialization.
        return typeof(value)
    end
    return value
end

function priming_configuration_key(record::ExecutionRecord, example::ExampleCase)
    options = Tuple(
        key => _canonical_priming_option(key, get(record.parameters, key, default)) for
        (key, default) in PRIMING_CODE_PATH_DEFAULTS
    )
    return (Int(example.variety["dim"]), options)
end

function _priming_configuration_label(key)
    dimension, options = key
    option_text = join(["$(name)=$(repr(value))" for (name, value) in options], ", ")
    return "(dimension=$(dimension), options={$(option_text)})"
end

function priming_record_errors(names::Vector{String}, records::Vector{ExecutionRecord})
    available = Set(record.name for record in records)
    unknown = sort!(unique([name for name in names if !(name in available)]))
    isempty(unknown) && return String[]
    return [
        "test/e2e/config.toml: priming.records contains unknown record(s): $(join(repr.(unknown), ", "))",
    ]
end

function priming_coverage_errors(
    catalogue::ExampleCatalogue,
    records::Vector{ExecutionRecord},
    config::CatalogueConfig,
)
    records_by_name = Dict(record.name => record for record in records)
    primed_configurations = Set{Any}()
    for name in config.priming_records
        record = get(records_by_name, name, nothing)
        record === nothing && continue
        example = get(catalogue.examples, record.example, nothing)
        example === nothing && continue
        push!(primed_configurations, priming_configuration_key(record, example))
    end

    configurations = Dict{Any,Vector{ExecutionRecord}}()
    for record in records
        example = get(catalogue.examples, record.example, nothing)
        example === nothing && continue
        key = priming_configuration_key(record, example)
        push!(get!(configurations, key, ExecutionRecord[]), record)
    end

    errors = String[]
    for key in sort!(collect(keys(configurations)); by = repr)
        cheapest =
            first(sort(configurations[key]; by = record -> (record.runtime, record.name)))
        cheapest.runtime < config.ci_seconds || continue
        key in primed_configurations && continue
        push!(
            errors,
            "$(_priming_configuration_label(key)): cheapest record $(repr(cheapest.name)) " *
            "has runtime $(cheapest.runtime)s (< ci_seconds=$(config.ci_seconds)); " *
            "add a priming record to test/e2e/config.toml",
        )
    end
    return errors
end

# Validates `catalogue` + `records` against `config` and the executor
# registry. Returns a Dict{String,Any} mapping each execution record's name
# to its resolved expected value (computed once here, reused at run time),
# or throws a CatalogueValidationError listing every problem found.
function validate_catalogue(
    catalogue::ExampleCatalogue,
    records::Vector{ExecutionRecord},
    config::CatalogueConfig,
)
    errors = String[]
    expected = Dict{String,Any}()

    for example in values(catalogue.examples)
        loc = "$(example.source_path) (example $(repr(example.id)))"
        local notes_values
        try
            notes_values = Common.parse_notes_invariants(example.notes)
        catch e
            push!(errors, "$(loc): $(sprint(showerror, e))")
            continue
        end
        result_primes = try
            Set(Common.result_prime(result) for result in example.results)
        catch e
            push!(
                errors,
                "$(loc): result entries must have valid prime values: $(sprint(showerror, e))",
            )
            continue
        end
        for (prime_string, invariant_values) in notes_values
            prime = parse(BigInt, prime_string)
            if !(prime in result_primes)
                push!(
                    errors,
                    "$(loc): derham-invariants notes tag has prime $(prime_string) with no result in the case",
                )
            end
            for invariant in keys(invariant_values)
                executor = get(INVARIANT_REGISTRY, invariant, nothing)
                if executor === nothing
                    push!(
                        errors,
                        "$(loc): derham-invariants notes tag names unregistered invariant $(repr(invariant))",
                    )
                elseif getproperty(executor, :EXPECTED_SOURCE) == :lpoly
                    push!(
                        errors,
                        "$(loc): invariant $(repr(invariant)) is derived from Lpoly and must not have a derham-invariants notes value",
                    )
                end
            end
        end
    end

    for record in records
        loc = "$(record.source_path) (record $(repr(record.name)))"

        if !isfinite(record.runtime) || record.runtime < 0
            push!(
                errors,
                "$(loc): runtime must be a finite, nonnegative number, got $(record.runtime)",
            )
        end

        for tag in record.tags
            tag in config.allowed_tags || push!(
                errors,
                "$(loc): unknown tag $(repr(tag)); allowed tags are $(sort(collect(config.allowed_tags)))",
            )
        end

        executor = get(EXECUTOR_REGISTRY, record.type, nothing)
        if executor === nothing
            push!(
                errors,
                "$(loc): unknown type $(repr(record.type)); known types are $(sort(collect(keys(EXECUTOR_REGISTRY))))",
            )
            continue
        end

        example = get(catalogue.examples, record.example, nothing)
        if example === nothing
            push!(errors, "$(loc): references unknown example $(repr(record.example))")
            continue
        end

        try
            executor.validate_parameters(record.name, record.parameters)
        catch e
            push!(errors, "$(loc): $(sprint(showerror, e))")
        end

        try
            expected[record.name] = executor.expected(example, record)
        catch e
            push!(errors, "$(loc): $(sprint(showerror, e))")
        end
    end

    append!(errors, priming_record_errors(config.priming_records, records))
    append!(errors, priming_coverage_errors(catalogue, records, config))

    isempty(errors) || throw(CatalogueValidationError(errors))
    return expected
end
