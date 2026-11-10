@testset "Degree-one functional equation sign" begin
    @test DeRham.sign_fe(2, 1, 3, [2], [3, 1]) == 1
    @test DeRham.sign_fe(2, 1, 3, [2], [6, 1]) == -1
    @test_throws ArgumentError DeRham.sign_fe(2, 1, 3, [1], [3, 1])

    @test_throws ArgumentError DeRham.sign_fe(2, 1, 2, [2], [2, 1])
    @test DeRham.sign_fe(2, 1, 2, [3], [2, 1]) == 1
    @test DeRham.sign_fe(2, 1, 2, [3], [6, 1]) == -1

    @test_throws ArgumentError DeRham.sign_fe(2, 1, 5, [2], [1, 1])
    @test_throws ArgumentError DeRham.sign_fe(2, 1, 5, [2], [10, 1])
end
