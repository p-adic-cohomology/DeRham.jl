#module PolynomialWithPole

#using Oscar

"""
gets all of the terms of polpol,
including their pole information

polpol - a PolynomialWithPole
"""
function wholeterms(polpol)
    t = Oscar.terms(polpol[1])
    result = []
    for i in t
        push!(result, [i, polpol[2]])
    end
    return result
end

"""
gets all of the terms of polpol,
forgetting their pole information
"""
numeratorterms(polpol) = Oscar.terms(polpol[1])


"""
Get the PolynomialWithPole which has the highest
pole order in the list

polpols - array of polynomials with pole
"""
function highestpoleorder(polpols)
    function cmp_poleorder(p1, p2)
        if p1[2] ≤ p2[2]
            p2
        else
            p1
        end
    end

    reduce(cmp_poleorder, polpols)[2]
end

"""
Gets all terms of a certain pole order in an
array of polynomials with pole

polpols - array of polynomials with pole
"""
function termsoforder(polpols, order)
    result = []
    for polpol in polpols
        if polpol[2] == order
            append!(result, wholeterms(polpol))
        end
    end

    result
end

#end
