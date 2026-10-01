@testset "newton_polygon" begin
    R, (x, y, z) = polynomial_ring(GF(7), ["x", "y", "z"])
    f = y^2 * z - x^3 - x * z^2 - z^3
    q = 7

    coeffs = DeRham.zeta_coefficients(f)

    @testset "data form matches the f form" begin
        @test DeRham.newton_polygon(coeffs, q) == DeRham.newton_polygon(f)
    end

    @testset "twisted-input round trip" begin
        p5_q4K3_sparse_fk_001 = [
            -476837158203125,
            95367431640625,
            15258789062500,
            -6103515625000,
            305175781250,
            61035156250,
            36621093750,
            -7324218750,
            -1220703125,
            439453125,
            0,
            0,
            -3515625,
            390625,
            93750,
            -18750,
            -1250,
            -250,
            200,
            -20,
            -5,
            1,
        ]
        q_k3 = 5

        for (coeffs, q) in [(coeffs, q), (p5_q4K3_sparse_fk_001, q_k3)]
            tw = DeRham.tate_twist(coeffs, q)
            normtw = DeRham.normalized_tate_twist(coeffs, q, 1)
            np = DeRham.newton_polygon(coeffs, q)
            @test DeRham.newton_polygon(tw, q) == np
            @test DeRham.newton_polygon(normtw, q) == np
        end
    end

    @testset "non-smooth returns false" begin
        R5, (y1, y2, y3, y4) = polynomial_ring(GF(5), ["x1", "x2", "x3", "x4"])
        f_nonsmooth = (y1^2 + y2^2)^2 + y3^4 + y4^4
        @test DeRham.newton_polygon(f_nonsmooth) == false
    end
end
