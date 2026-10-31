# test/e2e/test/runtests.jl — ordinary automated tests for the E2E harness
# itself (selection, capability handling, and validation logic). These run
# as part of the conventional test suite, independent of --workflow=/
# --name=/--tag= selection of the E2E catalogue's own cases, and do not
# execute any zeta_coefficients computation.

using Test

@testset "E2E harness" begin
    @testset "selection" begin
        @test E2E.parse_selection(String[]) == E2E.Selection(:workflow, :ci)
        @test E2E.parse_selection(["--workflow=presubmit"]) ==
              E2E.Selection(:workflow, :presubmit)
        @test E2E.parse_selection(["--workflow=full"]) == E2E.Selection(:workflow, :full)
        @test E2E.parse_selection(["--name=foo"]) == E2E.Selection(:name, "foo")
        @test E2E.parse_selection(["--tag=gpu"]) == E2E.Selection(:tag, "gpu")

        @test_throws ErrorException E2E.parse_selection(["--workflow=nonsense"])
        @test_throws ErrorException E2E.parse_selection(["--bogus=1"])
        @test_throws ErrorException E2E.parse_selection(["--name=a", "--tag=b"])
        @test_throws ErrorException E2E.parse_selection(["--workflow=ci", "--name=a"])

        cfg = E2E.CatalogueConfig(
            5.0,
            90.0,
            Set(["gpu", "multithreading", "uint", "bigint", "karatsuba"]),
        )
        mk(name, runtime; tags = String[]) = E2E.ExecutionRecord(
            name,
            "zeta_coefficients",
            tags,
            runtime,
            "ex",
            nothing,
            Dict{String,Any}(),
            "synthetic",
        )
        recs = [
            mk("fast", 1.0),
            mk("mid", 10.0),
            mk("slow", 200.0),
            mk("gpuslow", 200.0; tags = ["gpu"]),
        ]

        presubmit = E2E.select_records(recs, cfg, E2E.Selection(:workflow, :presubmit))
        ci = E2E.select_records(recs, cfg, E2E.Selection(:workflow, :ci))
        full = E2E.select_records(recs, cfg, E2E.Selection(:workflow, :full))
        presubmit_names = Set(r.name for r in presubmit.records)
        ci_names = Set(r.name for r in ci.records)
        full_names = Set(r.name for r in full.records)
        @test presubmit_names ⊆ ci_names
        @test ci_names ⊆ full_names
        @test presubmit_names == Set(["fast"])
        @test ci_names == Set(["fast", "mid"])
        @test full_names == Set(["fast", "mid", "slow", "gpuslow"])
        @test presubmit.broad && ci.broad && full.broad

        name_sel = E2E.select_records(recs, cfg, E2E.Selection(:name, "slow"))
        @test length(name_sel.records) == 1 && name_sel.records[1].name == "slow"
        @test !name_sel.broad
        @test_throws ErrorException E2E.select_records(
            recs,
            cfg,
            E2E.Selection(:name, "nope"),
        )

        tag_sel = E2E.select_records(recs, cfg, E2E.Selection(:tag, "gpu"))
        @test Set(r.name for r in tag_sel.records) == Set(["gpuslow"])
        @test tag_sel.broad
        @test_throws ErrorException E2E.select_records(
            recs,
            cfg,
            E2E.Selection(:tag, "not-a-tag"),
        )
        @test_throws ErrorException E2E.select_records(
            recs,
            cfg,
            E2E.Selection(:tag, "multithreading"),
        )
    end

    @testset "runtime summary" begin
        @test !E2E.runtime_drift(10.0, 5.0)
        @test !E2E.runtime_drift(10.0, 20.0)
        @test E2E.runtime_drift(10.0, 4.9)
        @test E2E.runtime_drift(10.0, 20.1)
        @test E2E.runtime_drift(0.0, 0.1)
        @test !E2E.runtime_drift(
            E2E.CaseOutcome("skip", false, true, "gpu required", 10.0, 0.0),
        )

        io = IOBuffer()
        E2E.print_runtime_summary(
            [
                E2E.CaseOutcome("slow", true, false, nothing, 1.0, 2.1),
                E2E.CaseOutcome("gpu_case", false, true, "gpu required", 1.0, 0.0),
            ],
            io,
        )
        summary = String(take!(io))
        @test occursin("declared (s)", summary)
        @test occursin("measured (s)", summary)
        @test occursin("slow", summary) && occursin("DRIFT", summary)
        @test occursin("gpu_case", summary) && occursin("SKIP", summary)
    end

    @testset "zeta_coefficients executor validation" begin
        exec = E2E.ZetaCoefficientsExecutor

        good_variety(
            kind = "plane_curve",
            vars = ["x1", "x2", "x3"];
            dim = 1,
            middle_factor_content = "full",
        ) = Dict{String,Any}(
            "coeff_domain" => Dict("kind" => "integer"),
            "dim" => dim,
            "non_middle_factors" => Dict(
                "kind" => "projective_lefschetz",
                "middle_factor_content" => middle_factor_content,
            ),
            "model" => Dict(
                "kind" => kind,
                "vars" => vars,
                "monomials" => [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                "coeffs" => [1, 1, 1],
            ),
        )
        example(variety, results) =
            E2E.ExampleCase("ex1", "synthetic.json", variety, results, "")
        record(; p = nothing, params = Dict{String,Any}()) = E2E.ExecutionRecord(
            "rec1",
            "zeta_coefficients",
            String[],
            0.1,
            "ex1",
            p,
            params,
            "synthetic.toml",
        )

        # correct ascending orientation round-trips to DeRham's descending
        # convention.
        ex = example(
            good_variety(),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
        )
        @test exec.expected(ex, record()) == BigInt[7, 1, 1]
        primitive_curve = example(
            good_variety(; middle_factor_content = "primitive"),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
        )
        @test exec.expected(primitive_curve, record()) == BigInt[7, 1, 1]

        full_threefold = example(
            good_variety(
                "projective_hypersurface",
                ["x1", "x2", "x3", "x4", "x5"];
                dim = 3,
            ),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
        )
        @test exec.expected(full_threefold, record()) == BigInt[7, 1, 1]

        primitive_surface = example(
            good_variety(
                "projective_hypersurface",
                ["x1", "x2", "x3", "x4"];
                dim = 2,
                middle_factor_content = "primitive",
            ),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 0, 0, 0, 0, 0, 7]))],
        )
        @test length(exec.expected(primitive_surface, record())) == 7

        full_surface = example(
            good_variety("projective_hypersurface", ["x1", "x2", "x3", "x4"]; dim = 2),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 0, 0, 0, 0, 0, 7]))],
        )
        @test_throws r"primitive middle factor" exec.expected(full_surface, record())

        # reversed (descending) coeffs_asc is caught as the specific
        # orientation mistake, not a generic error.
        bad = example(
            good_variety(),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [7, 1, 1]))],
        )
        @test_throws r"orientation mistake" exec.expected(bad, record())

        # unsupported model kind / coeff domain / non_middle_factors convention.
        unsupported_kind = example(
            good_variety("hyperelliptic"),
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
        )
        @test_throws ErrorException exec.expected(unsupported_kind, record())

        v = good_variety()
        v["coeff_domain"] = Dict("kind" => "number_field")
        unsupported_domain =
            example(v, [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))])
        @test_throws ErrorException exec.expected(unsupported_domain, record())

        v2 = good_variety()
        v2["non_middle_factors"] = Dict("kind" => "explicit")
        unsupported_nmf = example(v2, [Dict("p" => 7, "l_factors" => Dict())])
        @test_throws ErrorException exec.expected(unsupported_nmf, record())

        # multi-result example requires an explicit p on the record.
        multi = example(
            good_variety(),
            [
                Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7])),
                Dict("p" => 11, "Lpoly" => Dict("coeffs_asc" => [1, 0, 11])),
            ],
        )
        @test_throws ErrorException exec.expected(multi, record())
        @test exec.expected(multi, record(; p = 11)) == BigInt[11, 0, 1]

        # parameter validation.
        @test isnothing(
            exec.validate_parameters(
                "rec",
                Dict("algorithm" => "depthfirst", "S" => [0, 1]),
            ),
        )
        @test_throws ErrorException exec.validate_parameters(
            "rec",
            Dict("algorithm" => "not-an-algorithm"),
        )
        @test_throws ErrorException exec.validate_parameters(
            "rec",
            Dict("not_a_real_param" => 1),
        )

        # capability requirement derivation.
        @test :gpu in
              E2E.Common.required_capabilities(record(; params = Dict("use_gpu" => true)))
        @test isempty(E2E.Common.required_capabilities(record()))
        gpu_tagged = E2E.ExecutionRecord(
            "rec2",
            "zeta_coefficients",
            ["gpu"],
            0.1,
            "ex1",
            nothing,
            Dict{String,Any}(),
            "synthetic.toml",
        )
        @test :gpu in E2E.Common.required_capabilities(gpu_tagged)
        @test exec.INVARIANT_NAME == "zeta_coefficients"
        @test exec.EXPECTED_SOURCE == :lpoly
    end

    @testset "derived invariant executors" begin
        function variety(kind, dimension, vars, monomials, coefficients)
            return Dict{String,Any}(
                "coeff_domain" => Dict("kind" => "integer"),
                "dim" => dimension,
                "non_middle_factors" => Dict(
                    "kind" => "projective_lefschetz",
                    "middle_factor_content" => "primitive",
                ),
                "model" => Dict(
                    "kind" => kind,
                    "vars" => vars,
                    "monomials" => monomials,
                    "coeffs" => coefficients,
                ),
            )
        end
        function example_case(id, v, p, coeffs_desc)
            return E2E.ExampleCase(
                id,
                "synthetic.json",
                v,
                [Dict("p" => p, "Lpoly" => Dict("coeffs_asc" => reverse(coeffs_desc)))],
                "",
            )
        end
        function record(type; params = Dict{String,Any}())
            return E2E.ExecutionRecord(
                "record_$(type)",
                type,
                String[],
                0.1,
                "synthetic",
                nothing,
                params,
                "synthetic.toml",
            )
        end

        curve = example_case(
            "curve",
            variety(
                "plane_curve",
                1,
                ["x1", "x2", "x3"],
                [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                [1, 1, 1],
            ),
            7,
            BigInt[7, 1, 1],
        )
        curve_input_record = E2E.ExecutionRecord(
            "curve_input",
            "newton_polygon",
            String[],
            0.1,
            "curve",
            nothing,
            Dict{String,Any}(),
            "synthetic.toml",
        )
        curve_F = E2E.Common.build_input(curve, curve_input_record)
        newton = E2E.NewtonPolygonExecutor
        expected_newton = newton.expected(curve, curve_input_record)
        actual_newton = newton.normalize(newton.execute(curve_F, Dict{String,Any}()))
        @test newton.compare(expected_newton, actual_newton)
        @test !newton.compare(expected_newton, E2E.DeRham.SlopesPolygon([0, 2]))
        @test newton.INVARIANT_NAME == "newton_polygon"
        @test newton.EXPECTED_SOURCE == :lpoly
        @test isnothing(newton.validate_parameters("newton", Dict{String,Any}()))

        point_counts = E2E.PointCountsExecutor
        point_params = Dict{String,Any}("k" => 2)
        point_record = record("point_counts"; params = point_params)
        point_expected = point_counts.expected(curve, point_record)
        point_actual = point_counts.normalize(point_counts.execute(curve_F, point_params))
        @test point_counts.compare(point_expected, point_actual)
        @test point_expected == BigInt[9, 63]
        @test isnothing(point_counts.validate_parameters("counts", point_params))
        @test_throws r"requires parameter 'k'" point_counts.validate_parameters(
            "counts",
            Dict{String,Any}(),
        )
        @test_throws r"positive integer" point_counts.validate_parameters(
            "counts",
            Dict{String,Any}("k" => 0),
        )
        @test_throws r"unknown parameter" point_counts.validate_parameters(
            "counts",
            Dict{String,Any}("k" => 2, "other" => 1),
        )
        unsupported_curve = deepcopy(curve)
        unsupported_curve.variety["model"]["kind"] = "hyperelliptic"
        @test_throws r"does not support model.kind" newton.expected(
            unsupported_curve,
            curve_input_record,
        )
        @test_throws r"does not support model.kind" point_counts.expected(
            unsupported_curve,
            point_record,
        )

        cubic_surface = example_case(
            "cubic_surface",
            variety(
                "projective_hypersurface",
                2,
                ["x1", "x2", "x3", "x4"],
                [[3, 0, 0, 0], [0, 3, 0, 0], [0, 0, 3, 0], [0, 0, 0, 3]],
                [1, 1, 1, 1],
            ),
            11,
            BigInt[-1771561, 0, 43923, 0, -363, 0, 1],
        )
        bound = E2E.PicardRankBoundExecutor
        bound_record = record("picard_rank_bound")
        @test bound.expected(cubic_surface, bound_record) == 7
        @test bound.EXPECTED_SOURCE == :lpoly
        @test isnothing(bound.validate_parameters("bound", Dict{String,Any}()))
        @test_throws r"requires dimension 2" bound.expected(curve, bound_record)
        malformed_surface = deepcopy(cubic_surface)
        malformed_surface.variety["model"]["kind"] = "plane_curve"
        @test_throws r"requires model.kind .*projective_hypersurface" bound.expected(
            malformed_surface,
            bound_record,
        )

        k3_desc = vcat(BigInt[-7^21], zeros(BigInt, 20), BigInt[1])
        k3 = example_case(
            "k3",
            variety(
                "projective_hypersurface",
                2,
                ["x1", "x2", "x3", "x4"],
                [[4, 0, 0, 0], [0, 4, 0, 0], [0, 0, 4, 0], [0, 0, 0, 4]],
                [1, 1, 1, 1],
            ),
            7,
            k3_desc,
        )
        height = E2E.ArtinMazurHeightExecutor
        geo = E2E.K3GeometricPicardRankExecutor
        rank = E2E.K3PicardRankExecutor
        realization = E2E.K3PicardRealizationDegreeExecutor
        @test height.expected(k3, record("artin_mazur_height")) == Inf
        @test geo.expected(k3, record("k3_geometric_picard_rank")) == 22
        rank_params = Dict{String,Any}("k" => 2)
        @test rank.expected(k3, record("k3_picard_rank"; params = rank_params)) == 2
        @test isnothing(rank.validate_parameters("rank", rank_params))
        @test_throws r"requires parameter 'k'" rank.validate_parameters(
            "rank",
            Dict{String,Any}(),
        )
        @test_throws r"positive integer" rank.validate_parameters(
            "rank",
            Dict{String,Any}("k" => 0),
        )
        @test realization.expected(k3, record("k3_picard_realization_degree")) == 21
        @test_throws r"homogeneous degree-4" geo.expected(
            cubic_surface,
            record("k3_geometric_picard_rank"),
        )
        @test_throws r"homogeneous degree-4" height.expected(
            cubic_surface,
            record("artin_mazur_height"),
        )
        @test_throws r"homogeneous degree-4" rank.expected(
            cubic_surface,
            record("k3_picard_rank"; params = rank_params),
        )
        @test_throws r"homogeneous degree-4" realization.expected(
            cubic_surface,
            record("k3_picard_realization_degree"),
        )
        malformed_k3 = deepcopy(k3)
        malformed_k3.variety["model"]["kind"] = "plane_curve"
        @test_throws r"requires model.kind .*projective_hypersurface" geo.expected(
            malformed_k3,
            record("k3_geometric_picard_rank"),
        )
        for executor in (newton, bound, geo, realization, height)
            @test_throws r"unknown parameter" executor.validate_parameters(
                "unsupported",
                Dict{String,Any}("not_a_parameter" => 1),
            )
        end
        @test_throws r"unknown parameter" point_counts.validate_parameters(
            "unsupported",
            Dict{String,Any}("k" => 2, "not_a_parameter" => 1),
        )
        @test_throws r"unknown parameter" rank.validate_parameters(
            "unsupported",
            Dict{String,Any}("k" => 2, "not_a_parameter" => 1),
        )
        @test all(
            executor.EXPECTED_SOURCE == :lpoly for
            executor in (bound, geo, rank, realization, height)
        )
    end

    @testset "a_number notes executor" begin
        executor = E2E.ANumberExecutor
        @test executor.INVARIANT_NAME == "a_number"
        @test executor.EXPECTED_SOURCE == :notes
        @test E2E.INVARIANT_REGISTRY["a_number"] === executor
        @test isnothing(executor.validate_parameters("a_number_notes", Dict{String,Any}()))
        @test_throws r"unknown parameter" executor.validate_parameters(
            "a_number_notes",
            Dict{String,Any}("verbose" => 1),
        )

        variety = Dict{String,Any}(
            "coeff_domain" => Dict("kind" => "integer"),
            "dim" => 1,
            "non_middle_factors" => Dict(
                "kind" => "projective_lefschetz",
                "middle_factor_content" => "full",
            ),
            "model" => Dict(
                "kind" => "plane_curve",
                "vars" => ["x1", "x2", "x3"],
                "monomials" => [[0, 2, 1], [3, 0, 0], [1, 0, 2], [0, 0, 3]],
                "coeffs" => [1, -1, -1, -1],
            ),
        )
        function example(notes)
            return E2E.ExampleCase(
                "a_number_curve",
                "synthetic.json",
                deepcopy(variety),
                [Dict("p" => 7)],
                notes,
            )
        end
        record = E2E.ExecutionRecord(
            "a_number_notes",
            "a_number",
            String[],
            0.1,
            "a_number_curve",
            7,
            Dict{String,Any}(),
            "synthetic.toml",
        )

        valid = example("derham-invariants: {\"7\":{\"a_number\":0}}")
        @test executor.expected(valid, record) == 0
        cfg = E2E.CatalogueConfig(
            5.0,
            90.0,
            Set(["gpu", "multithreading", "uint", "bigint", "karatsuba"]),
        )
        expected = E2E.validate_catalogue(
            E2E.ExampleCatalogue(Dict(valid.id => valid)),
            [record],
            cfg,
        )
        @test expected[record.name] == 0

        @test_throws r"notes have no derham-invariants values" executor.expected(
            example("smooth plane curve"),
            record,
        )
        @test_throws r"malformed derham-invariants tag" executor.expected(
            example("derham-invariants: {broken}"),
            record,
        )

        curve = E2E.Common.build_input(valid, record)
        actual = executor.normalize(executor.execute(curve, record.parameters))
        @test executor.compare(expected[record.name], actual)
    end

    @testset "catalogue validation aggregates errors" begin
        cfg = E2E.CatalogueConfig(
            5.0,
            90.0,
            Set(["gpu", "multithreading", "uint", "bigint", "karatsuba"]),
        )
        variety = Dict{String,Any}(
            "coeff_domain" => Dict("kind" => "integer"),
            "dim" => 1,
            "non_middle_factors" => Dict(
                "kind" => "projective_lefschetz",
                "middle_factor_content" => "full",
            ),
            "model" => Dict(
                "kind" => "plane_curve",
                "vars" => ["x1", "x2", "x3"],
                "monomials" => [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                "coeffs" => [1, 1, 1],
            ),
        )
        good_example = E2E.ExampleCase(
            "ok",
            "s.json",
            variety,
            [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
            "",
        )
        catalogue = E2E.ExampleCatalogue(Dict("ok" => good_example))

        records = [
            E2E.ExecutionRecord(
                "r_ok",
                "zeta_coefficients",
                String[],
                0.1,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_missing_ref",
                "zeta_coefficients",
                String[],
                0.1,
                "does_not_exist",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_bad_tag",
                "zeta_coefficients",
                ["not-a-tag"],
                0.1,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_bad_type",
                "not-a-type",
                String[],
                0.1,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
            E2E.ExecutionRecord(
                "r_bad_runtime",
                "zeta_coefficients",
                String[],
                -1.0,
                "ok",
                nothing,
                Dict{String,Any}(),
                "s.toml",
            ),
        ]

        err = nothing
        try
            E2E.validate_catalogue(catalogue, records, cfg)
            error("expected CatalogueValidationError")
        catch e
            err = e
        end
        @test err isa E2E.CatalogueValidationError
        @test length(err.errors) == 4
        joined = join(err.errors, "\n")
        @test occursin("r_missing_ref", joined)
        @test occursin("r_bad_tag", joined)
        @test occursin("r_bad_type", joined)
        @test occursin("r_bad_runtime", joined)

        expected = E2E.validate_catalogue(catalogue, [records[1]], cfg)
        @test expected["r_ok"] == BigInt[7, 1, 1]

        validate_notes(notes) = begin
            ex = E2E.ExampleCase(
                "notes",
                "notes.json",
                deepcopy(variety),
                [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
                notes,
            )
            note_record = E2E.ExecutionRecord(
                "r_notes",
                "zeta_coefficients",
                String[],
                0.1,
                "notes",
                nothing,
                Dict{String,Any}(),
                "notes.toml",
            )
            E2E.validate_catalogue(
                E2E.ExampleCatalogue(Dict("notes" => ex)),
                [note_record],
                cfg,
            )
        end

        valid_notes = "prose\nderham-invariants: {\"7\":{\"a_number\":0}}\n"
        parsed_notes = E2E.Common.parse_notes_invariants(valid_notes)
        @test parsed_notes == Dict("7" => Dict("a_number" => 0))
        @test_throws r"canonical prime string" E2E.Common.parse_notes_invariants(
            "derham-invariants: {\"9\":{\"a_number\":1}}",
        )
        @test E2E.Common.notes_invariant_value(
            E2E.ExampleCase(
                "notes",
                "notes.json",
                deepcopy(variety),
                [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
                valid_notes,
            ),
            7,
            "a_number",
        ) == 0
        @test_throws r"malformed derham-invariants tag" validate_notes(
            "derham-invariants: {broken}",
        )
        @test_throws r"malformed derham-invariants tag" validate_notes(
            "derham-invariants {\"7\":{}}",
        )
        @test_throws r"only one tag line is allowed" validate_notes(
            "derham-invariants: {}\nderham-invariants: {}",
        )
        @test_throws r"prime 11 with no result in the case" validate_notes(
            "derham-invariants: {\"11\":{\"a_number\":1}}",
        )
        @test_throws r"unregistered invariant" validate_notes(
            "derham-invariants: {\"7\":{\"not_registered\":1}}",
        )
        @test_throws r"derived from Lpoly" validate_notes(
            "derham-invariants: {\"7\":{\"zeta_coefficients\":[7,1,1]}}",
        )
    end

    @testset "artin_invariant notes executor" begin
        executor = E2E.ArtinInvariantExecutor
        @test executor.INVARIANT_NAME == "artin_invariant"
        @test executor.EXPECTED_SOURCE == :notes
        @test E2E.INVARIANT_REGISTRY["artin_invariant"] === executor
        @test isnothing(
            executor.validate_parameters("artin_invariant_notes", Dict{String,Any}()),
        )
        @test_throws r"artin_invariant accepts no parameters" executor.validate_parameters(
            "artin_invariant_notes",
            Dict{String,Any}("S" => [0, 1, 2]),
        )

        examples_path = joinpath(
            @__DIR__,
            "..",
            "examples",
            "zenodo",
            "zenodo_22784396_quartic_k3.json",
        )
        cases = E2E._read_example_document(examples_path)
        parse_case(id) = E2E._parse_example_case(
            only(filter(case -> case["id"] == id, cases)),
            examples_path,
            1,
        )
        with_notes(example, notes) = E2E.ExampleCase(
            example.id,
            example.source_path,
            example.variety,
            example.results,
            notes,
        )
        sigma_case = parse_case("zenodo_22784396_p5_r6")
        ordinary_case = parse_case("zenodo_22784396_p5_r172")
        function artin_record(example; p = nothing, params = Dict{String,Any}())
            return E2E.ExecutionRecord(
                "artin_invariant_$(example.id)",
                "artin_invariant",
                String[],
                0.1,
                example.id,
                p,
                params,
                "synthetic.toml",
            )
        end

        record = artin_record(sigma_case)
        @test executor.expected(sigma_case, record) == 9
        @test executor.compare(9, executor.normalize(9))
        @test_throws r"requires an odd prime field" executor.expected(
            sigma_case,
            artin_record(sigma_case; p = 2),
        )
        @test_throws r"requires a supersingular quartic K3" executor.expected(
            ordinary_case,
            artin_record(ordinary_case; p = 5),
        )

        missing_tag =
            with_notes(sigma_case, "Zenodo 22784396 quartic K3 without its invariant note")
        @test_throws r"no derham-invariants values" executor.expected(missing_tag, record)
        malformed_tag = with_notes(sigma_case, "derham-invariants: {broken}")
        @test_throws r"malformed derham-invariants tag" executor.expected(
            malformed_tag,
            record,
        )
    end

    @testset "loader — example documents" begin
        good_variety = Dict{String,Any}(
            "coeff_domain" => Dict("kind" => "integer"),
            "dim" => 1,
            "non_middle_factors" => Dict(
                "kind" => "projective_lefschetz",
                "middle_factor_content" => "full",
            ),
            "model" => Dict(
                "kind" => "plane_curve",
                "vars" => ["x1", "x2", "x3"],
                "monomials" => [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
                "coeffs" => [1, 1, 1],
            ),
        )
        good_case(id) = Dict{String,Any}(
            "id" => id,
            "variety" => good_variety,
            "results" => [Dict("p" => 7, "Lpoly" => Dict("coeffs_asc" => [1, 1, 7]))],
            "notes" => "",
        )

        mktempdir() do dir
            bad = joinpath(dir, "bad.json")
            write(bad, "{ not valid json")
            @test_throws r"failed to parse as JSON" E2E._read_example_document(bad)

            notobj = joinpath(dir, "notobj.json")
            write(notobj, "[1, 2, 3]")
            @test_throws r"top level must be a JSON object" E2E._read_example_document(
                notobj,
            )

            wrongver = joinpath(dir, "wrongver.json")
            write(wrongver, """{"schema_version": "2", "cases": []}""")
            @test_throws r"schema_version must be" E2E._read_example_document(wrongver)

            nocases = joinpath(dir, "nocases.json")
            write(nocases, """{"schema_version": "3"}""")
            @test_throws r"missing required key 'cases'" E2E._read_example_document(nocases)

            casesnotarr = joinpath(dir, "casesnotarr.json")
            write(casesnotarr, """{"schema_version": "3", "cases": {}}""")
            @test_throws r"'cases' must be an array" E2E._read_example_document(casesnotarr)

            good = joinpath(dir, "good.json")
            write(good, """{"schema_version": "3", "cases": [{"id": "x"}]}""")
            cases = E2E._read_example_document(good)
            @test length(cases) == 1
            @test cases[1]["id"] == "x"
        end

        @test E2E._parse_example_case(good_case("ex1"), "p.json", 1).id == "ex1"

        missing_key = Dict{String,Any}(
            "id" => "ex1",
            "variety" => good_variety,
            "results" => [Dict("p" => 7)],
        )
        @test_throws r"missing required key 'notes'" E2E._parse_example_case(
            missing_key,
            "p.json",
            1,
        )

        empty_id = merge(good_case("ex1"), Dict{String,Any}("id" => ""))
        @test_throws r"empty or non-string id" E2E._parse_example_case(
            empty_id,
            "p.json",
            1,
        )

        bad_variety_type = merge(good_case("ex1"), Dict{String,Any}("variety" => "nope"))
        @test_throws r"'variety' must be an object" E2E._parse_example_case(
            bad_variety_type,
            "p.json",
            1,
        )

        missing_variety_key = merge(
            good_case("ex1"),
            Dict{String,Any}("variety" => Dict{String,Any}("coeff_domain" => Dict())),
        )
        @test_throws r"variety is missing required key" E2E._parse_example_case(
            missing_variety_key,
            "p.json",
            1,
        )

        empty_results = merge(good_case("ex1"), Dict{String,Any}("results" => []))
        @test_throws r"'results' must be a nonempty array" E2E._parse_example_case(
            empty_results,
            "p.json",
            1,
        )

        non_array_results =
            merge(good_case("ex1"), Dict{String,Any}("results" => "notarray"))
        @test_throws r"'results' must be a nonempty array" E2E._parse_example_case(
            non_array_results,
            "p.json",
            1,
        )

        case_json(id) = """
        {"id": "$(id)", "variety": {"coeff_domain": {"kind": "integer"}, "dim": 1, "non_middle_factors": {"kind": "projective_lefschetz", "middle_factor_content": "full"}, "model": {"kind": "plane_curve", "vars": ["x1", "x2", "x3"], "monomials": [[3, 0, 0], [0, 3, 0], [0, 0, 3]], "coeffs": [1, 1, 1]}}, "results": [{"p": 7, "Lpoly": {"coeffs_asc": [1, 1, 7]}}], "notes": ""}
        """
        document_json(cases_json) = """{"schema_version": "3", "cases": [$(cases_json)]}"""

        mktempdir() do dir
            dup_within = joinpath(dir, "dup_within.json")
            write(dup_within, document_json(case_json("dup") * ", " * case_json("dup")))
            into = Dict{String,E2E.ExampleCase}()
            @test_throws r"duplicate example id" E2E._load_examples_into!(into, dup_within)

            f1 = joinpath(dir, "a.json")
            write(f1, document_json(case_json("shared")))
            f2 = joinpath(dir, "b.json")
            write(f2, document_json(case_json("shared")))
            into2 = Dict{String,E2E.ExampleCase}()
            E2E._load_examples_into!(into2, f1)
            @test_throws r"duplicate example id" E2E._load_examples_into!(into2, f2)
        end
    end

    @testset "loader — execution records" begin
        mktempdir() do dir
            records_dir = joinpath(dir, "records")
            mkpath(records_dir)

            bad = joinpath(records_dir, "bad.toml")
            write(bad, "name = \"unterminated")
            @test_throws r"failed to parse as TOML" E2E.load_execution_records(dir)
            rm(bad)

            missing_key = joinpath(records_dir, "missing_key.toml")
            write(missing_key, "[[tests]]\nname = \"x\"\n")
            @test_throws r"missing required key 'test'" E2E.load_execution_records(dir)
            rm(missing_key)

            not_array = joinpath(records_dir, "not_array.toml")
            write(not_array, "test = \"nope\"\n")
            @test_throws r"'test' must be an array of tables" E2E.load_execution_records(
                dir,
            )
            rm(not_array)

            missing_field = joinpath(records_dir, "missing_field.toml")
            write(missing_field, "[[test]]\nname = \"r1\"\ntype = \"zeta_coefficients\"\n")
            @test_throws r"missing required key 'example'" E2E.load_execution_records(dir)
            rm(missing_field)

            empty_name = joinpath(records_dir, "empty_name.toml")
            write(
                empty_name,
                "[[test]]\nname = \"\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\n",
            )
            @test_throws r"empty or non-string name" E2E.load_execution_records(dir)
            rm(empty_name)

            r1 = joinpath(records_dir, "r1.toml")
            write(
                r1,
                "[[test]]\nname = \"dup\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\n",
            )
            r2 = joinpath(records_dir, "r2.toml")
            write(
                r2,
                "[[test]]\nname = \"dup\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\n",
            )
            @test_throws r"duplicate execution record name" E2E.load_execution_records(dir)
            rm(r1)
            rm(r2)

            ok = joinpath(records_dir, "ok.toml")
            write(
                ok,
                "[[test]]\nname = \"ok1\"\ntype = \"zeta_coefficients\"\nexample = \"ex\"\nruntime = 1.0\n",
            )
            recs = E2E.load_execution_records(dir)
            @test length(recs) == 1 && recs[1].name == "ok1"
        end
    end
end
