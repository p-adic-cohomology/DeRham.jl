using DeRham, Oscar, Test, SHA
function compositions(total, parts)
    parts==1 && return [[total]]
    return [vcat(a, tail) for a = 0:total for tail in compositions(total-a, parts-1)]
end
function main(path)
    counts=Dict{Int,Int}()
    refusals=0
    for line in eachline(path)
        fields=split(line, '\t'; keepempty = true)
        n, p, seed=parse.(Int, fields[1:3])
        status=fields[4]
        d=n==4 ? 4 : 3
        coefficients=parse.(Int, split(fields[5], ','))
        R, xs=polynomial_ring(GF(p), ["x$i" for i = 1:n])
        f=sum(R(c)*prod(xs .^ e) for (e, c) in zip(compositions(d, n), coefficients))
        if status=="ok"
            expected=parse.(BigInt, split(fields[6], ','))
            @test BigInt.(reverse(DeRham.dense_dwork_coefficients(f)))==expected
            counts[n]=get(counts, n, 0)+1
        elseif status=="singular"
            @test_throws ArgumentError DeRham.dense_dwork_coefficients(f)
            refusals+=1
        else
            @test_throws ErrorException DeRham.dense_dwork_coefficients(f)
            refusals+=1
        end
    end
    println(
        "cohort matches=",
        sort(collect(counts)),
        " refusals=",
        refusals,
        " sha256=",
        bytes2hex(sha256(read(path))),
    )
end
main(only(ARGS))
