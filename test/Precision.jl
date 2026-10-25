
function _reference_ibm_floor_log(p, a)
    value = BigInt(a)
    base = BigInt(p)
    exponent = 0
    while value >= base
        value = div(value, base)
        exponent += 1
    end
    exponent
end

function _reference_ibm_f0(p, m, n)
    sum(_reference_ibm_floor_log(p, max(big(1), BigInt(m) - i)) for i in 1:n)
end

function _reference_vp(value, p)
    remaining = BigInt(value)
    prime = BigInt(p)
    exponent = 0
    while rem(remaining, prime) == 0
        remaining = div(remaining, prime)
        exponent += 1
    end
    exponent
end

function _reference_ibm_nn(p, m, n, bounds, K, l)
    index = (BigInt(m) + l) * p
    A = index <= K ? bounds[Int(index)] : _reference_ibm_f0(p, index, n)
    g = _reference_vp(binomial(BigInt(m) + l - 1, l), p)
    A - l - g
end

function _reference_ibm_constant(p, m, n, bounds, K)
    best = _reference_ibm_nn(p, m, n, bounds, K, big(1))
    L = big(2)
    while true
        exponent = L + best
        if BigInt(m) + L >= 2n + 1 &&
           exponent >= 0 &&
           (BigInt(p) * (BigInt(m) + L) - 1)^n <= BigInt(p)^exponent
            return best
        end
        best = max(best, _reference_ibm_nn(p, m, n, bounds, K, L))
        L += 1
    end
end

function _reference_ibm_crank_bound_h(p, n, K)
    bounds = BigInt[_reference_ibm_f0(p, k, n) for k in 1:K]
    for m in 1:div(K, p)
        pm = p * m
        N = _reference_ibm_constant(p, m, n, bounds, K)
        updated = min(bounds[pm], max(n - 1 + bounds[m], N))
        bounds[pm] = updated
        for k in (pm - p + 1):pm
            bounds[k] = min(bounds[k], updated)
        end
    end
    bounds
end

function test_akr_lattice_bound()
    for (p, k) in ((3, 5), (3, 10), (11, 7))
        power = big(p)^k
        @test DeRham._ibm_floor_log(p, power) == k
        @test DeRham._ibm_floor_log(p, power - 1) == k - 1
    end

    @test DeRham.ibm_derham_bound(3, 244, 2) == 9
    @test DeRham.ibm_lattice_correction(2, 3) == 1
    @test DeRham.ibm_lattice_correction(2, 5) == 3
    @test DeRham.ibm_lattice_correction(5, 3) == 0

    @test DeRham._ibm_crank_bound_h(2, 2, 4) == [0, 0, 1, 2]
    @test DeRham.ibm_crank_bound(2, 2, 4) == [0, 0, 1, 2]
    @test DeRham._ibm_crank_bound_h(2, 3, 6) == [0, 0, 1, 2, 4, 5]
    @test DeRham.ibm_crank_bound(2, 3, 6) == [1, 1, 2, 3, 5, 6]

    tail_bounds = [DeRham.ibm_derham_bound(2, k, 1) for k in 1:120]
    tail_bounds[116] = 2
    tail_bounds[118] = 3
    tail_bounds[120] = 3
    tail_values = [
        (2 * (57 + l) <= 120 ? tail_bounds[2 * (57 + l)] :
         DeRham.ibm_derham_bound(2, 2 * (57 + l), 1)) - l -
        _reference_vp(binomial(big(57 + l - 1), big(l)), 2) for l in 1:6
    ]
    @test tail_values == BigInt[1, 1, 0, 2, 1, 0]
    @test DeRham.ibm_constant(2, 57, tail_bounds, 1) == 2

    regression_cases = [
        (2, 3, 6),
        (2, 10, 1000),
        (3, 2, 732),
        (3, 2, 12),
        (2, 1, 120),
        (5, 2, 11),
        (101, 3, 303),
    ]
    for (p, n, K) in regression_cases
        H = DeRham._ibm_crank_bound_h(p, n, K)
        f0 = [DeRham.ibm_derham_bound(p, k, n) for k in 1:K]
        reference = _reference_ibm_crank_bound_h(p, n, K)
        c = DeRham.ibm_lattice_correction(p, n)

        @test length(H) == K
        @test all(H .<= f0)
        @test BigInt.(H) == reference
        @test BigInt.(DeRham.ibm_crank_bound(p, n, K)) == reference .+ c
    end

    for p in (2, 3, 5, 101), n in (2, 3), multiplier in (1, 2, 3)
        K = p * n * multiplier
        H = DeRham._ibm_crank_bound_h(p, n, K)
        f0 = [DeRham.ibm_derham_bound(p, k, n) for k in 1:K]
        @test all(H .<= f0)
        @test BigInt.(H) == _reference_ibm_crank_bound_h(p, n, K)
    end
end

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
