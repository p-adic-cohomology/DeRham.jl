@testset "pchunk threaded execution" begin
    if Threads.nthreads() < 2
        if get(ENV, "DERHAM_PCHUNK_THREAD_TEST_CHILD", "0") == "1"
            @test Threads.nthreads() >= 2
        else
            test_dir = normpath(joinpath(@__DIR__, ".."))
            e2e_module = joinpath(test_dir, "e2e", "E2E.jl")
            threading_test = abspath(@__FILE__)
            test_environment = dirname(Base.active_project())
            child_script = """
                using Test
                using DeRham
                include($(repr(e2e_module)))
                using .E2E
                include($(repr(threading_test)))
            """
            command = `$(Base.julia_cmd()) --threads=2 --project=$test_environment -e $child_script`
            run(addenv(command, "DERHAM_PCHUNK_THREAD_TEST_CHILD" => "1"))
        end
    else
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
end
