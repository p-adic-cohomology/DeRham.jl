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
    boat_p5_sparse = DeRham.boat_shape_Lpoly(p5_q4K3_sparse_fk_001, 21, 5)
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
    boat_fermat_p13 = DeRham.boat_shape_Lpoly(fermat_quartic_p13, 21, 13)
    @test DeRham.artin_mazur_height(fermat_quartic_p13, 13) == 1
    @test DeRham.artin_mazur_height(boat_fermat_p13, 13) == 1
end
