@testset "Picard rank" begin
    @testset "_is_boat_shape" begin
        @test DeRham._is_boat_shape([1, 2, 5], 5) == true
        @test DeRham._is_boat_shape([1, 2, 7], 5) == false
        @test DeRham._is_boat_shape([5, 2, 1], 5) == false
        @test DeRham._is_boat_shape([-5, 2, 5], 5) == true
    end

    @testset "_cyclotomic_orders" begin
        fermat_cubic_p11 = [-1771561, 0, 43923, 0, -363, 0, 1]
        @test DeRham._cyclotomic_orders(fermat_cubic_p11, 11) == Dict(1 => 3, 2 => 3)

        p5_q4K3_dense_001 = [
            -476837158203125,
            -114440917968750,
            7629394531250,
            7629394531250,
            1220703125000,
            -152587890625,
            -109863281250,
            -15869140625,
            2197265625,
            1123046875,
            117187500,
            -23437500,
            -8984375,
            -703125,
            203125,
            56250,
            3125,
            -1000,
            -250,
            -10,
            6,
            1,
        ]
        @test DeRham._cyclotomic_orders(p5_q4K3_dense_001, 5) == Dict(1 => 1)

        boat_shape_p5_dense = DeRham.boat_shape_Lpoly(p5_q4K3_dense_001, 21, 5)
        @test DeRham._cyclotomic_orders(boat_shape_p5_dense, 5) == Dict(1 => 1)
    end

    @testset "coefficient API" begin
        q_generic, p5_q4K3_dense_001 = 5,
        [
            -476837158203125,
            -114440917968750,
            7629394531250,
            7629394531250,
            1220703125000,
            -152587890625,
            -109863281250,
            -15869140625,
            2197265625,
            1123046875,
            117187500,
            -23437500,
            -8984375,
            -703125,
            203125,
            56250,
            3125,
            -1000,
            -250,
            -10,
            6,
            1,
        ]
        @test DeRham.k3_geometric_picard_rank(p5_q4K3_dense_001, q_generic) == 2
        @test DeRham.k3_picard_realization_degree(p5_q4K3_dense_001, q_generic) == 1
        @test DeRham.k3_picard_rank(p5_q4K3_dense_001, q_generic, 1) == 2

        q_row3, k3f3_row3 = 3,
        [
            10460353203,
            -4649045868,
            0,
            0,
            86093442,
            28697814,
            -19131876,
            1594323,
            -1062882,
            708588,
            -118098,
            -39366,
            26244,
            -4374,
            729,
            -972,
            162,
            54,
            0,
            0,
            -4,
            1,
        ]
        @test DeRham.k3_geometric_picard_rank(k3f3_row3, q_row3) == 4
        @test DeRham.k3_picard_realization_degree(k3f3_row3, q_row3) == 2
        @test DeRham.k3_picard_rank(k3f3_row3, q_row3, 1) == 3
        @test DeRham.k3_picard_rank(k3f3_row3, q_row3, 2) == 4

        q_row5, k3f3_row5 = 3,
        [
            -10460353203,
            2324522934,
            -387420489,
            -129140163,
            43046721,
            -14348907,
            0,
            1594323,
            -1062882,
            0,
            118098,
            -39366,
            0,
            4374,
            -729,
            0,
            81,
            -27,
            9,
            3,
            -2,
            1,
        ]
        @test DeRham.k3_geometric_picard_rank(k3f3_row5, q_row5) == 8
        @test DeRham.k3_picard_realization_degree(k3f3_row5, q_row5) == 7
        @test DeRham.k3_picard_rank(k3f3_row5, q_row5, 1) == 2
        @test DeRham.k3_picard_rank(k3f3_row5, q_row5, 7) == 8

        q_sparse, p5_q4K3_sparse_fk_001 = 5,
        [
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
        @test DeRham.k3_geometric_picard_rank(p5_q4K3_sparse_fk_001, q_sparse) == 14
        @test DeRham.k3_picard_realization_degree(p5_q4K3_sparse_fk_001, q_sparse) == 24
        for (k, rk) in
            [1 => 4, 2 => 6, 3 => 6, 4 => 6, 6 => 10, 8 => 10, 12 => 10, 24 => 14]
            @test DeRham.k3_picard_rank(p5_q4K3_sparse_fk_001, q_sparse, k) == rk
        end

        boat_sparse = DeRham.boat_shape_Lpoly(p5_q4K3_sparse_fk_001, 21, q_sparse)
        @test boat_sparse[1] == -5
        @test DeRham.k3_geometric_picard_rank(boat_sparse, q_sparse) == 14
        @test DeRham.k3_picard_realization_degree(boat_sparse, q_sparse) == 24

        q_fermat4, fermat_quartic_p13 = 13,
        [
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
        @test DeRham.k3_geometric_picard_rank(fermat_quartic_p13, q_fermat4) == 20
        @test DeRham.k3_picard_realization_degree(fermat_quartic_p13, q_fermat4) == 2
        @test DeRham.k3_picard_rank(fermat_quartic_p13, q_fermat4, 1) == 8
        @test DeRham.k3_picard_rank(fermat_quartic_p13, q_fermat4, 2) == 20

        q_row2, k3f3_row2 = 3,
        [
            -10460353203,
            0,
            1162261467,
            -387420489,
            -258280326,
            43046721,
            28697814,
            -9565938,
            -4782969,
            1062882,
            531441,
            -177147,
            -39366,
            19683,
            4374,
            -1458,
            -243,
            162,
            27,
            -9,
            0,
            1,
        ]
        @test DeRham.k3_geometric_picard_rank(k3f3_row2, q_row2) == 22
        @test DeRham.k3_picard_realization_degree(k3f3_row2, q_row2) == 12
        for (k, rk) in [1 => 2, 2 => 4, 3 => 6, 4 => 4, 6 => 14, 12 => 22]
            @test DeRham.k3_picard_rank(k3f3_row2, q_row2, k) == rk
        end

        fermat_cubic_p11 = [-1771561, 0, 43923, 0, -363, 0, 1]
        @test DeRham.picard_rank_bound(fermat_cubic_p11, 11) == 7

        fermat_quintic_p7 = [
            -88124787089723195184393736687912818113311201,
            0,
            0,
            0,
            477143786824823630735992743416437582454413,
            0,
            0,
            0,
            -1192362649291520943113684490003592459278,
            0,
            0,
            0,
            1820906447619982003922050449818064286,
            0,
            0,
            0,
            -1895987554789652232322001717844715,
            0,
            0,
            0,
            1421398416751925871794920071687,
            0,
            0,
            0,
            -789336896710773773036190516,
            0,
            0,
            0,
            328753393049051967112116,
            0,
            0,
            0,
            -102692646725026645287,
            0,
            0,
            0,
            23761545357264715,
            0,
            0,
            0,
            -3958608139486,
            0,
            0,
            0,
            449654478,
            0,
            0,
            0,
            -31213,
            0,
            0,
            0,
            1,
        ]
        @test DeRham.picard_rank_bound(fermat_quintic_p7, 7) == 53

        @test_throws ArgumentError DeRham.k3_geometric_picard_rank(fermat_cubic_p11, 11)
        @test_throws ArgumentError DeRham.k3_picard_realization_degree(fermat_cubic_p11, 11)
        @test_throws ArgumentError DeRham.k3_picard_rank(fermat_cubic_p11, 11, 1)
    end

    @testset "polynomial wrappers" begin
        R3, (x1, x2, x3, x4) = polynomial_ring(GF(3), ["x1", "x2", "x3", "x4"])
        f_fermat3 = x1^4 + x2^4 + x3^4 + x4^4

        @test DeRham.picard_rank_bound(f_fermat3) == 22
        @test DeRham.k3_geometric_picard_rank(f_fermat3) == 22
        @test DeRham.k3_picard_realization_degree(f_fermat3) == 2
        @test DeRham.k3_picard_rank(f_fermat3, 1) == 12
        @test DeRham.k3_picard_rank(f_fermat3, 2) == 22

        R5, (y1, y2, y3, y4) = polynomial_ring(GF(5), ["x1", "x2", "x3", "x4"])
        f_nonsmooth = (y1^2 + y2^2)^2 + y3^4 + y4^4
        @test DeRham.zeta_coefficients(f_nonsmooth) == false
        @test DeRham.picard_rank_bound(f_nonsmooth) == false
        @test DeRham.k3_geometric_picard_rank(f_nonsmooth) == false

        @test_throws ArgumentError DeRham.k3_geometric_picard_rank(
            y1^3 + y2^3 + y3^3 + y4^3,
        )
    end
end
