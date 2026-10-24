using DeRham, Oscar, Pkg, SHA, TOML

function rational_coordinate(ev)
    exponents = copy(ev)
    pole = 5
    coordinate = big(1)//big(1)
    while any(!iszero, exponents)
        index = findfirst(!iszero, exponents)
        exponent = exponents[index]
        exponent == 1 && return big(0)//big(1)
        coordinate *= (big(exponent)-1)//(2big(pole-1))
        exponents[index] -= 2
        pole -= 1
    end
    @assert pole == 3
    return coordinate
end

function reproduce_computeT(output_path)
    field = GF(3)
    polynomial_ring_fp, variables = polynomial_ring(field, ["x$i" for i in 0:5])
    f = sum(x^2 for x in variables)
    g = variables[1]^2 * variables[2]^2
    params = DeRham.default_params()
    cache = DeRham.PolyExpCache(6, params.termorder, 4)
    bases = DeRham.compute_monomial_bases(f, params, cache)
    @assert length.(bases) == [0,0,1,0,0]
    @assert bases[3] == [one(polynomial_ring_fp)]
    monomials = DeRham.gen_mon(cache[4], polynomial_ring_fp)
    column = findfirst(==(g), monomials)
    @assert cache[4][column] == [2,2,0,0,0,0]
    exact_coordinates = rational_coordinate.(cache[4])
    expected = factorial(big(4)) .* exact_coordinates
    @assert exact_coordinates[column] == big(1)//big(48)
    @assert expected[column] == big(1)//big(2)
    source_method = which(DeRham.computeT, (typeof(f),typeof(bases),Int,typeof(params),typeof(cache)))
    probes = []
    for precision in [2,3,4]
        actual = DeRham.computeT(f,bases,precision,params,cache)
        residue = base_ring(actual)
        expected_mod = [residue(numerator(c))*inv(residue(denominator(c))) for c in expected]
        unweighted_mod = [residue(numerator(c))*inv(residue(denominator(c))) for c in expected .// 2]
        @assert size(actual) == (1,126)
        @assert actual[1,column] == inv(residue(4))
        @assert actual[1,column] != inv(residue(2))
        @assert all(actual[1,j] == unweighted_mod[j] for j in eachindex(expected))
        @assert all(2actual[1,j] == expected_mod[j] for j in eachindex(expected))
        push!(probes, Dict(
            "M"=>precision,"modulus"=>string(modulus(residue)),
            "shape"=>collect(size(actual)),"g_column_one_based"=>column,
            "observed_g"=>string(actual[1,column]),"expected_g"=>string(expected_mod[column]),
            "observed_equals_one_quarter"=>true,"expected_equals_one_half"=>true,
            "all_columns_equal_half_expected"=>true,"all_columns_factor_two_recovers_expected"=>true,
            "mismatching_columns"=>count(j->actual[1,j]!=expected_mod[j],eachindex(expected)),
            "actual_row"=>[string(actual[1,j]) for j in eachindex(expected)],
            "expected_row"=>string.(expected_mod)))
    end
    dependency_info = []
    for (uuid,info) in sort!(collect(Pkg.dependencies());by=x->x[2].name)
        push!(dependency_info,Dict("name"=>info.name,"uuid"=>string(uuid),
            "version"=>string(info.version),"tree_hash"=>string(info.tree_hash),
            "source"=>string(info.source),"is_tracking_path"=>info.is_tracking_path))
    end
    result = Dict(
        "assignment"=>"dr-eevm.9","VERSION"=>string(VERSION),"word_size"=>Sys.WORD_SIZE,
        "project"=>Base.active_project(),"derham_path"=>pathof(DeRham),"oscar_path"=>pathof(Oscar),
        "image"=>unsafe_string(Base.JLOptions().image_file),"image_id"=>get(ENV,"KAIMON_CITY_SYSIMAGE","unset"),
        "pid"=>getpid(),"n"=>5,"d"=>2,"p"=>3,"F"=>string(f),"g"=>string(g),
        "basis_lengths"=>length.(bases),"basis_pole3"=>string.(bases[3]),"termorder"=>string(params.termorder),
        "computeT_method_file"=>string(source_method.file),"computeT_method_line"=>source_method.line,
        "exponent_vectors"=>cache[4],
        "exact_g_coordinate"=>string(exact_coordinates[column]),"exact_Tg"=>string(expected[column]),
        "source_sha256"=>bytes2hex(sha256(read(joinpath(dirname(pathof(DeRham)),"FinalReduction.jl")))),
        "project_sha256"=>bytes2hex(sha256(read(Base.active_project()))),
        "manifest_sha256"=>bytes2hex(sha256(read(joinpath(dirname(Base.active_project()),"Manifest-v1.13.toml")))),
        "probes"=>probes,"dependency_info"=>dependency_info)
    open(output_path,"w") do io
        TOML.print(io,result)
    end
    return (version=VERSION,loaded_source=pathof(DeRham),basis_lengths=length.(bases),
        coordinate=exact_coordinates[column],expected_T=expected[column],
        probes=[(M=p["M"],modulus=p["modulus"],observed=p["observed_g"],expected=p["expected_g"],mismatches=p["mismatching_columns"]) for p in probes],
        output=output_path)
end

reproduce_computeT(joinpath(@__DIR__,"native-results.toml"))
