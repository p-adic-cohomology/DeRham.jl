function akr_relation_identities_hold(table, f)
    partials = [derivative(f, i) for i = 1:3]
    for (monomial, cofactors) in table.relations
        reconstructed = Dict{NTuple{3,Int},BigInt}()
        for (i, cofactor, relation_coefficient) in cofactors
            for term in terms(partials[i])
                exponent = ntuple(j -> cofactor[j] + exponent_vector(term, 1)[j], Val(3))
                coefficient = mod(
                    get(reconstructed, exponent, big(0)) +
                    relation_coefficient * DeRham._akr_bigint(lift(ZZ, coeff(term, 1))),
                    table.modulus,
                )
                if iszero(coefficient)
                    delete!(reconstructed, exponent)
                else
                    reconstructed[exponent] = coefficient
                end
            end
        end
        reconstructed == Dict(monomial => big(1)) || return false
    end
    return true
end

@testset "AKR cutoff relation table" begin
    _, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    f = x^4 + y^4 + z^4 + x*y*z*(x+y+z)
    table = DeRham._build_akr_relation_table(f, 7, 1)

    @test table.cutoff == 7
    @test length(table.relations) == 36
    @test akr_relation_identities_hold(table, f)
end

@testset "AKR cutoff BigInt reduction map" begin
    precision_modulus = big(7)^23
    precision_ring, = residue_ring(ZZ, precision_modulus)
    _, (x, y, z) = polynomial_ring(precision_ring, ["x", "y", "z"])
    f = x^4 + y^4 + z^4 + x*y*z*(x+y+z)
    table = DeRham._build_akr_relation_table(f, 7, 23)
    @test table.modulus > typemax(UInt64)
    @test akr_relation_identities_hold(table, f)

    m = 3
    params = (; termorder = :invlex)
    source_degree = 4 * m - 3
    cofactor_degree = source_degree - 4 + 1
    source_monomials = DeRham.gen_exp_vec(3, source_degree, params.termorder)
    cofactor_monomials = DeRham.gen_exp_vec(3, cofactor_degree, params.termorder)
    cofactor_indices = Dict(Tuple(e) => i for (i, e) in enumerate(cofactor_monomials))
    g = [precision_ring(i) for i in eachindex(source_monomials)]
    expected = fill(precision_ring(0), 3 * length(cofactor_monomials))

    for (j, exponent) in enumerate(source_monomials)
        iszero(g[j]) && continue
        left = table.cutoff
        h = ntuple(Val(3)) do i
            part = min(left, exponent[i])
            left -= part
            part
        end
        shift = ntuple(i -> exponent[i] - h[i], Val(3))
        for (i, cofactor, factor) in table.relations[h]
            shifted = ntuple(l -> shift[l] + cofactor[l], Val(3))
            position = (i - 1) * length(cofactor_monomials) + cofactor_indices[shifted]
            expected[position] += g[j] * precision_ring(factor)
        end
    end

    actual = DeRham._akr_cofactor_vector(g, m, f, table, params)
    @test actual == expected
    cofactor_length = length(cofactor_monomials)
    expected_reduced_polynomial = derivative(
        DeRham.vector_to_polynomial(
            expected[1:cofactor_length],
            2,
            cofactor_degree,
            parent(f),
            params.termorder,
        ),
        1,
    )
    for i = 1:2
        component = DeRham.vector_to_polynomial(
            expected[(i*cofactor_length+1):((i+1)*cofactor_length)],
            2,
            cofactor_degree,
            parent(f),
            params.termorder,
        )
        expected_reduced_polynomial += derivative(component, i + 1)
    end
    expected_reduced =
        DeRham.polynomial_to_vector(expected_reduced_polynomial, 3, params.termorder)
    @test DeRham._reduce_akr_level(g, m, f, table, params) == expected_reduced

    cached_map = DeRham._akr_bigint_reduction_map(table, m, f, params.termorder)
    @test cached_map === DeRham._akr_bigint_reduction_map(table, m, f, params.termorder)
    @test length(table.bigint_reduction_maps) == 1
end

@testset "AKR UInt64 reduction path" begin
    p = 7
    precision = 12
    modulus = big(p)^precision
    precision_ring, = residue_ring(ZZ, modulus)
    _, (x, y, z) = polynomial_ring(precision_ring, ["x", "y", "z"])
    f = x^3 + y^3 + z^3
    table = DeRham._build_akr_relation_table(f, p, precision)

    @test table.modulus <= typemax(UInt64)

    m = 3
    params = (; termorder = :invlex)
    source_monomials = DeRham.gen_exp_vec(3, 3 * m - 3, params.termorder)
    g = fill(precision_ring(modulus - 1), length(source_monomials))
    reduction_map = DeRham._akr_reduction_map(table, m, f, params.termorder)
    source_residue = UInt64(modulus - 1)
    @test any(
        UInt128(source_residue) * UInt128(reduction_map.factors[k]) > typemax(UInt64) for
        j in eachindex(g) for k = reduction_map.offsets[j]:(reduction_map.offsets[j+1]-1)
    )

    cofactor_degree = 3 * m - 3 - 3 + 1
    cofactor_monomials = DeRham.gen_exp_vec(3, cofactor_degree, params.termorder)
    cofactor_indices = Dict(Tuple(e) => i for (i, e) in enumerate(cofactor_monomials))
    expected_cofactors_bigint = fill(big(0), 3 * length(cofactor_monomials))
    for (j, exponent) in enumerate(source_monomials)
        left = table.cutoff
        h = ntuple(Val(3)) do i
            part = min(left, exponent[i])
            left -= part
            part
        end
        @test left == 0
        shift = ntuple(i -> exponent[i] - h[i], Val(3))
        for (i, cofactor, factor) in table.relations[h]
            shifted = ntuple(l -> shift[l] + cofactor[l], Val(3))
            position = (i - 1) * length(cofactor_monomials) + cofactor_indices[shifted]
            expected_cofactors_bigint[position] =
                mod(expected_cofactors_bigint[position] + (modulus - 1) * factor, modulus)
        end
    end

    actual_cofactors = DeRham._akr_cofactor_vector(g, m, f, table, params)
    @test [DeRham._akr_bigint(lift(ZZ, c)) for c in actual_cofactors] == expected_cofactors_bigint

    expected_cofactors = precision_ring.(expected_cofactors_bigint)
    expected_reduced =
        DeRham._akr_pseudoinverse_derivative(expected_cofactors, m, f, [0, 1, 2], params)
    @test DeRham._reduce_akr_level(g, m, f, table, params) == expected_reduced
end

@testset "AKR partial-variable fallback" begin
    p = 2
    coefficient_ring, = residue_ring(ZZ, p)
    _, (x, y, z) = polynomial_ring(coefficient_ring, ["x", "y", "z"])
    f = x^3 + y^3 + z^3
    n = 2
    d = 3
    S = [0, 1]
    params = DeRham.default_params()
    cache = DeRham.controlled_reduction_cache(n, d, S, params)
    l = d * n - n + d - length(S)
    pseudo_inverse = DeRham.pseudo_inverse_controlled_lifted(f, S, l, 1, params, cache)
    picache = DeRham._build_akr_picache(n, n, d, S, f, pseudo_inverse, cache, params, p)

    @test !(picache isa DeRham._AKRRelationTable)
    @test picache[n] === pseudo_inverse

    curve = DeRham.fermat_hypersurface(3, 3, p)
    expected = BigInt[2, 0, 1]
    for algorithm in (:akr_naive, :akr_levelbylevel)
        @test DeRham.zeta_coefficients(curve; algorithm = algorithm, S = S) == expected
    end
end

@testset "AKR full-variable fallback when characteristic divides degree" begin
    p = 2
    coefficient_ring, = residue_ring(ZZ, p)
    _, (x, y, z) = polynomial_ring(coefficient_ring, ["x", "y", "z"])
    f = x^3 * y + y^3 * z + z^3 * x
    n = 2
    d = 4
    S = [0, 1, 2]
    params = DeRham.default_params()
    pseudo_inverse = Matrix{Int}(undef, 0, 0)
    picache = DeRham._build_akr_picache(n, n, d, S, f, pseudo_inverse, nothing, params, p)

    @test !(picache isa DeRham._AKRRelationTable)
    @test picache[n] === pseudo_inverse
end

@testset "AKR BigInt multi-level reduction" begin
    curve = DeRham.fermat_hypersurface(3, 3, 7)
    hodge_polygon, relative_precision, series_precision, _ =
        DeRham.calculate_precision_information(curve, DeRham.get_basis_of_cohomology(curve))
    big_precision = (hodge_polygon, relative_precision, series_precision, 24)
    @test big(7)^24 > typemax(UInt64)

    expected = BigInt[7, 1, 1]
    for algorithm in (:akr_naive, :akr_levelbylevel)
        @test DeRham.zeta_coefficients(
            curve;
            algorithm = algorithm,
            precision_info = big_precision,
        ) == expected
    end
end

@testset "AKR cutoff selectors on a Fermat cubic" begin
    cubic = DeRham.fermat_hypersurface(3, 3, 7)
    expected = BigInt[7, 1, 1]
    for algorithm in (:akr, :akr_naive, :akr_levelbylevel)
        @test DeRham.zeta_coefficients(cubic; algorithm = algorithm) == expected
    end
end
