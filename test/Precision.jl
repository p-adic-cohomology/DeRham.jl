
function test_hodge_polygon_values()
    k3 = DeRham.SlopesPolygon([1, 20, 1])

    @test k3[1] == 0
    @test k3[2] == 1
    @test k3[10] == 9
    @test k3[21] == 20
    @test k3[22] == 22

    genus5 = DeRham.SlopesPolygon([5, 5])

    @test genus5[0] == 0
    @test genus5[2] == 0
    @test genus5[5] == 0
    @test genus5[6] == 1
    @test genus5[8] == 3
    @test genus5[10] == 5
end

function test_hodge_polygon_examples()
    R, (x, y, z, w) = polynomial_ring(GF(7), 4)

    f = x^4 + y^4 + z^4 + w^4

    @test DeRham.hodge_polygon(f) == DeRham.SlopesPolygon([1, 19, 1])
end

function test_hodge_polygon_data_form_matches_f_form()
    p = 11
    for (n, d) in [(2, 3), (3, 3), (2, 4), (3, 4)]
        f = DeRham.fermat_hypersurface(n + 1, d, p)
        @test DeRham.hodge_polygon(n, d) == DeRham.hodge_polygon(f)
    end
end

function test_user_relative_precision_bypasses_default_table()
    p = 3
    _, (x, y, z, w) = polynomial_ring(GF(p), 4)
    f = x^4 + y^4 + z^4 + w^4
    basis = DeRham.get_basis_of_cohomology(f)
    r_m = [0, 20, 20]

    calculated = DeRham.calculate_series_precision(p, 3, r_m)
    expected_user = max.(calculated, 1)
    default_table = DeRham.series_precision_default(p, 3, 4)

    user_info = DeRham.precision_information_user_r_m(f, r_m; basis = basis)
    @test calculated[1] == 0
    @test user_info[3] == expected_user
    @test user_info[3] != default_table

    uniform_r = fill(20, 3)
    uniform_info = DeRham.precision_information_user_r(f, 20; basis = basis)
    expected_uniform = max.(DeRham.calculate_series_precision(p, 3, uniform_r), 1)
    @test uniform_info[3] == expected_uniform

    auto_info = DeRham.precision_information_auto(f; basis = basis)
    @test auto_info[3] == default_table

    legacy_auto_info = DeRham.precision_information(f, basis)
    @test legacy_auto_info[3] == default_table

    small_user = zeros(Int, 3)
    small_max_auto_info =
        DeRham.precision_information_max_auto_r_m(f, small_user; basis = basis)
    @test small_max_auto_info[2] == auto_info[2]
    @test small_max_auto_info[3] == auto_info[3]
    @test small_max_auto_info[4] == auto_info[4]

    max_auto_info = DeRham.precision_information_max_auto_r_m(f, r_m; basis = basis)
    r_m_max = max.(auto_info[2], r_m)
    calculated_max = DeRham.calculate_series_precision(p, 3, r_m_max)
    expected_max = max.(default_table, calculated_max)
    expected_max[expected_max .== 0] .= 1
    @test max_auto_info[2] == r_m_max
    @test max_auto_info[3] == expected_max
    @test all(max_auto_info[3] .>= calculated_max)
end

function test_max_auto_precision_without_default_table()
    p = 5
    _, (x, y, z, w) = polynomial_ring(GF(p), 4)
    f = x^3 + y^3 + z^3 + w^3
    basis = DeRham.get_basis_of_cohomology(f)

    @test DeRham.series_precision_default(p, 3, 3) === nothing

    auto_info = DeRham.precision_information_auto(f; basis = basis)
    user_r = zeros(Int, length(auto_info[2]))
    max_auto_info = DeRham.precision_information_max_auto_r_m(f, user_r; basis = basis)
    expected = max.(DeRham.calculate_series_precision(p, 3, auto_info[2]), 1)

    @test max_auto_info[2] == auto_info[2]
    @test max_auto_info[3] == expected
end
