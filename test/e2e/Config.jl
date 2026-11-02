# Config.jl — config.toml model: workflow runtime thresholds and the
# allowed-tag set. Tags have no implied relationships; directory placement
# has no execution semantics.

struct CatalogueConfig
    presubmit_seconds::Float64
    ci_seconds::Float64
    allowed_tags::Set{String}
    priming_records::Vector{String}
end

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

    haskey(raw, "tags") || error("$(path): missing [tags] table")
    haskey(raw["tags"], "allowed") || error("$(path): tags.allowed is required")
    allowed_tags = Set{String}(String.(raw["tags"]["allowed"]))

    haskey(raw, "priming") || error("$(path): missing [priming] table")
    haskey(raw["priming"], "records") || error("$(path): priming.records is required")
    raw_priming_records = raw["priming"]["records"]
    raw_priming_records isa AbstractVector ||
        error("$(path): priming.records must be an array of record names")
    all(name -> name isa AbstractString, raw_priming_records) ||
        error("$(path): priming.records must contain only strings")
    priming_records = String.(raw_priming_records)
    any(isempty, priming_records) &&
        error("$(path): priming.records must not contain empty names")
    length(unique(priming_records)) == length(priming_records) ||
        error("$(path): priming.records must not contain duplicate names")

    return CatalogueConfig(presubmit_seconds, ci_seconds, allowed_tags, priming_records)
end
