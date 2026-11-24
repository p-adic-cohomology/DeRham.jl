@testset "AKR cutoff relation table" begin
    _, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    f = x^4 + y^4 + z^4 + x*y*z*(x+y+z)
    table = DeRham._build_akr_relation_table(f, 7, 1)
    partials = [derivative(f, i) for i = 1:3]

    @test table.cutoff == 7
    @test length(table.relations) == 36
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
        @test reconstructed == Dict(monomial => big(1))
    end
end

@testset "AKR cutoff selectors on a Fermat cubic" begin
    cubic = DeRham.fermat_hypersurface(3, 3, 7)
    expected = BigInt[7, 1, 1]
    for algorithm in (:akr, :akr_naive, :akr_levelbylevel)
        @test DeRham.zeta_coefficients(cubic; algorithm = algorithm) == expected
    end
end
