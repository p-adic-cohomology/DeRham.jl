module DeRham

using Oscar:
    Oscar,
    GF,
    MatElem,
    MPolyBuildCtx,
    MPolyRingElem,
    QQ,
    QQFieldElem,
    ZZ,
    ZZMatrix,
    ZZModMatrix,
    ZZRingElem,
    base_ring,
    basis,
    change_base_ring,
    characteristic,
    charpoly,
    coeff,
    coefficient_ring,
    coefficients,
    data,
    derivative,
    divexact,
    euler_phi,
    evaluate,
    exponent_vector,
    factor,
    finish,
    fraction_field,
    gens,
    grade,
    hom,
    homogenizer,
    ideal,
    is_cyclotomic_polynomial_with_data,
    is_homogeneous,
    is_prime,
    is_smooth,
    leading_coefficient,
    lift,
    map_entries,
    matrix,
    matrix_space,
    modulus,
    ncols,
    nrows,
    number_of_rows,
    nvars,
    polynomial_ring,
    proj,
    push_term!,
    quo,
    residue_ring,
    rref,
    special_linear_group,
    terms,
    total_degree,
    valuation,
    vars,
    zero!,
    zero_matrix,
    zzModMatrix
using BitIntegers: BitIntegers
using LinearAlgebra: LinearAlgebra, I, nullspace, rank
using Combinatorics: Combinatorics
using Memoize: Memoize
using LRUCache: LRUCache, LRU
using LFUDACache: LFUDACache, LFUDA
using OhMyThreads: OhMyThreads

using CUDA: CUDA, @cuda, blockDim, blockIdx, threadIdx, totalmem
using GPUFiniteFieldMatrices:
    GPUFiniteFieldMatrices,
    CuModArray,
    CuModMatrix,
    CuModVector,
    KaratsubaArray,
    KaratsubaMatrix,
    KaratsubaVector
# Pkg.add(url="https://github.com/UCSD-computational-number-theory/GPUFiniteFieldMatrices.jl")

# comment this out when not debugging
#using Debugger

#include("NemoAdditions.jl")

include("Utils.jl")
include("ZetaFunctionParams.jl")
include("presentations/CohomologyPresentation.jl")
include("frontends/ZetaFrontend.jl")
include("frontends/ProjectiveHypersurface.jl")
include("graded_ring/GradedExpCache.jl")
include("recurrence/backend/NemoDerivedWrappers.jl")
include("recurrence/backend/LinearAlgebraWrappers.jl")
include("graded_ring/GradedPieces.jl")
include("graded_ring/PolynomialVectors.jl")
include("graded_ring/LinearAlgebra.jl")
include("graded_ring/maps/MonomialBasis.jl")
include("graded_ring/maps/ControlledMap.jl")
include("graded_ring/maps/SmoothNondegenerate.jl")
include("graded_ring/maps/StandardReduction.jl")
include("SlopesPolygon.jl")
include("reduction/PolynomialWithPole.jl")
include("precision/SeriesPrecisionDefaults.jl")
include("precision/PrecisionEstimate.jl")
include("reduction/ControlledReductionCache.jl")
include("recurrence/PEP.jl")
include("recurrence/EvaluatePEP.jl")
include("reduction/ComputeRuv.jl")
include("reduction/Context.jl")
include("reduction/Steps.jl")
include("recurrence/backend/GPU.jl")
include("recurrence/RuvStorage.jl")
include("reduction/policies/PChunk.jl")
include("reduction/policies/VarByVar.jl")
include("reduction/policies/DepthFirst.jl")
include("reduction/policies/AKR.jl")
include("reduction/ControlledReduction.jl")

include("frobenius/Frobenius.jl")
include("finalreduction/FinalReduction.jl")
include("CharPolyFrob.jl")

include("ExamplePolynomials.jl")

include("precision/PrecisionInformation.jl")
include("K3OverF3.jl")

include("./invariants/ZetaCoefficients.jl")
include("./invariants/HodgePolygon.jl")
include("./invariants/Smoothness.jl")
include("./invariants/NewtonPolygon.jl")
include("./invariants/HasseWittMatrix.jl")
include("./invariants/aNumber.jl")
include("./invariants/FrobeniusMatrix.jl")
include("./invariants/CohomologyBasis.jl")
include("./invariants/LPolynomial.jl")
include("./invariants/TateTwist.jl")
include("./invariants/PointCounts.jl")
include("./invariants/ZetaFunction.jl")
include("./invariants/PicardRank.jl")
include("./invariants/ArtinMazurHeight.jl")
include("./invariants/ArtinInvariant.jl")
include("./invariants/Jacobian.jl")


export zeta_coefficients,
    zeta_coefficients_with_precision,
    frobenius_matrix_with_precision,
    zeta_function,
    L_polynomial,
    point_counts,
    newton_polygon,
    hodge_polygon,
    hasse_witt_matrix,
    a_number,
    picard_rank_bound,
    k3_geometric_picard_rank,
    k3_picard_rank,
    k3_picard_realization_degree,
    artin_mazur_height,
    artin_invariant,
    jacobian_zeta_coefficients,
    jacobian_point_counts,
    jacobian_zeta_function

public fermat_hypersurface,
    random_hypersurface,
    check_smoothness,
    SlopesPolygon,
    tate_twist,
    tate_untwist,
    normalized_tate_twist,
    ibm_lattice_correction,
    pregen_default_context

end
