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
end
