@testset "degree-one L-polynomial sign" begin
    @test DeRham.sign_fe(2, 1, 3, BigInt[3, 2], BigInt[6, 1]; degree_one_sign_factor = 2) == 1
    @test DeRham.sign_fe(2, 1, 3, BigInt[3, 2], BigInt[-6, 1]; degree_one_sign_factor = 2) == -1
    @test DeRham.sign_fe(2, 1, 3, BigInt[2, 1], BigInt[3, 1]) == 1
    @test DeRham.sign_fe(2, 1, 3, BigInt[2, 1], BigInt[-3, 1]) == -1
    @test_throws ArgumentError DeRham.sign_fe(2, 1, 3, BigInt[1, 1], BigInt[3, 1])
    @test_throws ArgumentError DeRham.sign_fe(2, 1, 2, BigInt[2, 1], BigInt[4, 1]; degree_one_sign_factor = 2)

    @test DeRham._cyclotomic_orders(BigInt[3, 1], 3) == Dict(2 => 1)
    @test DeRham._cyclotomic_orders(BigInt[-3, 1], 3) == Dict(1 => 1)
end
