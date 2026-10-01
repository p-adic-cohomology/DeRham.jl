@testset "Artin-Mazur height" begin
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
    boat_p5_sparse = DeRham.normalized_tate_twist(p5_q4K3_sparse_fk_001, 5, 1)
    @test boat_p5_sparse[1] == -5
    @test DeRham.artin_mazur_height(p5_q4K3_sparse_fk_001, 5) == 2
    @test DeRham.artin_mazur_height(boat_p5_sparse, 5) == 2

    fermat_quartic_p13 = [
        -247064529073450392704413,
        -80405615970649536087235,
        -224910813903914786258,
        2508620616620588000570,
        188312900398839895003,
        -29841375627214911331,
        -3496390230500968632,
        172033060544399704,
        27724721295752390,
        -827978102045094,
        -111941095381388,
        8610853490876,
        376867593102,
        -74670735230,
        -2741627512,
        329708184,
        16651063,
        -621751,
        -49010,
        26,
        55,
        1,
    ]
    boat_fermat_p13 = DeRham.normalized_tate_twist(fermat_quartic_p13, 13, 1)
    @test DeRham.artin_mazur_height(fermat_quartic_p13, 13) == 1
    @test DeRham.artin_mazur_height(boat_fermat_p13, 13) == 1

    @testset "twisted-input round trip" begin
        q = 5
        tw = DeRham.tate_twist(p5_q4K3_sparse_fk_001, q)
        @test DeRham.artin_mazur_height(tw, q) ==
              DeRham.artin_mazur_height(p5_q4K3_sparse_fk_001, q)
    end

    @testset "f form" begin
        R5, (a, b, c, d) = polynomial_ring(GF(5), ["a", "b", "c", "d"])
        fk3 = a^4 + b^4 + c^4 + d^4
        @test DeRham.artin_mazur_height(fk3) ==
              DeRham.artin_mazur_height(DeRham.zeta_coefficients(fk3), 5)

        R5b, (y1, y2, y3, y4) = polynomial_ring(GF(5), ["x1", "x2", "x3", "x4"])
        f_nonsmooth = (y1^2 + y2^2)^2 + y3^4 + y4^4
        @test DeRham.artin_mazur_height(f_nonsmooth) == false
    end
end
