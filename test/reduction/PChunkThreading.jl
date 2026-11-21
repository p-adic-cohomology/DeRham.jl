@testset "pchunk threaded execution" begin
    e2e_root = joinpath(@__DIR__, "..", "e2e")
    catalogue = E2E.load_example_catalogue(e2e_root)
    records = E2E.load_execution_records(e2e_root)
    record = only(filter(item -> item.name == "pchunk__dim_1_deg_3_fermat_1", records))
    example = catalogue.examples[record.example]
    input = E2E.Common.build_input(example, record)
    parameters = Dict{String,Any}(record.parameters)
    executor = E2E.ZetaCoefficientsExecutor

    parameters["use_threads"] = false
    single_threaded = executor.execute(input, parameters)
    parameters["use_threads"] = true
    threaded = executor.execute(input, parameters)

    @test threaded == single_threaded
end
