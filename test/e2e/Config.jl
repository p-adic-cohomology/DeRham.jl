# Config.jl — config.toml model: workflow runtime thresholds and the
# allowed-tag set. Tags have no implied relationships; directory placement
# has no execution semantics.

struct CatalogueConfig
    presubmit_seconds::Float64
    ci_seconds::Float64
    allowed_tags::Set{String}
    priming_records::Vector{String}
    presubmit_budget_seconds::Float64
    presubmit_records::Union{Nothing,Vector{String}}
end

CatalogueConfig(presubmit_seconds, ci_seconds, allowed_tags, priming_records) =
    CatalogueConfig(presubmit_seconds, ci_seconds, allowed_tags, priming_records, Inf)

CatalogueConfig(presubmit_seconds, ci_seconds, allowed_tags, priming_records, budget) =
    CatalogueConfig(
        presubmit_seconds,
        ci_seconds,
        allowed_tags,
        priming_records,
        budget,
        nothing,
    )

is_presubmit_record(record::ExecutionRecord, config::CatalogueConfig) =
    record.runtime < config.presubmit_seconds &&
    (config.presubmit_records === nothing || record.name in config.presubmit_records)

function load_config(path::AbstractString)
    isfile(path) || error("E2E config not found: $(path)")
    raw = TOML.parsefile(path)

    haskey(raw, "thresholds") || error("$(path): missing [thresholds] table")
    thresholds = raw["thresholds"]
    haskey(thresholds, "presubmit_seconds") ||
        error("$(path): thresholds.presubmit_seconds is required")
    haskey(thresholds, "ci_seconds") || error("$(path): thresholds.ci_seconds is required")
    presubmit_seconds = Float64(thresholds["presubmit_seconds"])
    ci_seconds = Float64(thresholds["ci_seconds"])
    presubmit_seconds > 0 || error("$(path): thresholds.presubmit_seconds must be > 0")
    ci_seconds > presubmit_seconds ||
        error("$(path): thresholds.ci_seconds must be > thresholds.presubmit_seconds")
    presubmit_budget_seconds = Float64(get(thresholds, "presubmit_budget_seconds", Inf))
    presubmit_budget_seconds > 0 ||
        error("$(path): thresholds.presubmit_budget_seconds must be > 0")

    presubmit_records = nothing
    if haskey(raw, "presubmit")
        smoke = raw["presubmit"]
        smoke isa AbstractDict || error("$(path): [presubmit] must be a table")
        haskey(smoke, "records") || error("$(path): presubmit.records is required")
        names = smoke["records"]
        names isa AbstractVector ||
            error("$(path): presubmit.records must be an array of record names")
        all(name -> name isa AbstractString, names) ||
            error("$(path): presubmit.records must contain only strings")
        presubmit_records = String.(names)
        isempty(presubmit_records) && error("$(path): presubmit.records must not be empty")
        any(isempty, presubmit_records) &&
            error("$(path): presubmit.records must not contain empty names")
        length(unique(presubmit_records)) == length(presubmit_records) ||
            error("$(path): presubmit.records must not contain duplicate names")
    end

    haskey(raw, "tags") || error("$(path): missing [tags] table")
    haskey(raw["tags"], "allowed") || error("$(path): tags.allowed is required")
    allowed_tags = Set{String}(String.(raw["tags"]["allowed"]))

    haskey(raw, "priming") || error("$(path): missing [priming] table")
    priming = raw["priming"]
    priming isa AbstractDict || error("$(path): [priming] must be a table")
    haskey(priming, "records") || error("$(path): priming.records is required")
    raw_priming_records = priming["records"]
    raw_priming_records isa AbstractVector ||
        error("$(path): priming.records must be an array of record names")
    all(name -> name isa AbstractString, raw_priming_records) ||
        error("$(path): priming.records must contain only strings")
    priming_records = String.(raw_priming_records)
    any(isempty, priming_records) &&
        error("$(path): priming.records must not contain empty names")
    length(unique(priming_records)) == length(priming_records) ||
        error("$(path): priming.records must not contain duplicate names")

    return CatalogueConfig(
        presubmit_seconds,
        ci_seconds,
        allowed_tags,
        priming_records,
        presubmit_budget_seconds,
        presubmit_records,
    )
end
