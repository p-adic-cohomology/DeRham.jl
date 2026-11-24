@testset "pairs of quadrics" begin
    count_pair = DeRham.quadric_pencil_point_count
    zeta_pair = DeRham.quadric_pencil_zeta
    P, T = polynomial_ring(ZZ, "T")

    # On the conic [s^2:s*u:u^2] the second quadric cuts out
    # s^4 + s*u^3 + 2*u^4, irreducible over F_3. No point lies at infinity.
    K = GF(3)
    R, (x, y, z) = polynomial_ring(K, ["x", "y", "z"])
    f, g = x * z - y^2, x^2 + y * z + 2 * z^2
    @test count_pair(f, g) == 0
    @test zeta_pair(f, g) == one(P) // (1 - T^4)

    L, a = Oscar.finite_field(3, 2, "a")
    _, (u, v, w) = polynomial_ring(L, ["u", "v", "w"])
    f2, g2 = u * w - v^2, u^2 + v * w + 2 * w^2
    @test count_pair(f2, g2 + a * f2) == 0
    @test zeta_pair(f2, g2 + a * f2) == one(P) // (1 - T^2)^2

    # P^4 over F_7. Direct enumeration gives 64 common zeros. All five
    # singular pencil members are rational, hence three + signs and two -.
    _, X = polynomial_ring(GF(7), ["x$i" for i = 1:5])
    A = sum(X .^ 2)
    B = sum((i - 1) * X[i]^2 for i = 1:5)
    expected = one(P) // ((1 - T) * (1 - 7T)^4 * (1 + 7T)^2 * (1 - 49T))
    @test count_pair(A, B) == 64
    @test zeta_pair(A, B) == expected
    @test zeta_pair(B, A) == expected  # simple root in the infinity chart
    # Swap coordinates by a shear and change the generators of the pencil.
    Y = [X[1] + X[2]; X[2:end]]
    C, D = sum(Y .^ 2), sum((i - 1) * Y[i]^2 for i = 1:5)
    @test zeta_pair(D, C + D) == expected
    @test count_pair(C, D) == 64

    # In seven variables the discriminant is t^7-t and all seven signs
    # are +1 (each signed cofactor equals the negative derivative, +1).
    _, W = polynomial_ring(GF(7), ["w$i" for i = 1:7])
    ambient = prod(1 - ZZ(7)^i * T for i = 0:4)
    @test zeta_pair(sum(W .^ 2), sum((i - 1) * W[i]^2 for i = 1:7)) ==
          one(P) // (ambient * (1 - 49T)^7)

    # The sole common point is [0:0:1]. Interleaved cross terms exercise
    # rings whose display order differs from the internal term order.
    m = 2x^2 + x*z + y*z
    h = 2x^2 + x*y + 2x*z + y^2 + y*z
    @test count_pair(m, h) == 1

    @test count_pair(x^2, x^2) == 4
    @test count_pair(x * y, x * z) == 5  # a line and an isolated point
    @test_throws ArgumentError zeta_pair(x^2, x^2)
    @test_throws ArgumentError zeta_pair(x^2 + y^2 + z^2, x^2)
    @test_throws ArgumentError zeta_pair(x^2, x^2 + y^2 + z^2)
    @test_throws ArgumentError count_pair(zero(R), f)
    @test_throws ArgumentError count_pair(f + x, g)
    @test_throws ArgumentError count_pair(f, f2)
    _, U = polynomial_ring(GF(2), ["u", "v", "w"])
    @test_throws ArgumentError count_pair(U[1]^2, U[2]^2)
    _, V = polynomial_ring(GF(3), ["v$i" for i = 1:4])
    @test count_pair(V[1]^2, V[2]^2) == 4  # underlying projective line
    @test_throws ArgumentError zeta_pair(V[1]^2, V[2]^2)
end
