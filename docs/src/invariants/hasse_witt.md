# Hasse–Witt matrix & a-number

For a smooth plane curve `C` of genus `g` over `GF(p)`, the Hasse–Witt
matrix represents Frobenius on `H^1(C, O_C)` in the basis used by the
Frobenius computation. DeRham obtains it by reducing the final `g x g` block
of a precision-1 Frobenius matrix modulo `p`. See the
[Frobenius matrix page](frobenius_matrix.md) for that matrix and its
precision.

The curve's a-number is `g - rank(H)`, where `H` is the Hasse–Witt matrix.
For an elliptic curve, a full-rank matrix gives the ordinary case and a zero
matrix gives the supersingular case.

The Fermat cubic over `GF(7)` is ordinary because 7 is congruent to 1 modulo
3. The Fermat cubic over `GF(5)` is supersingular because 5 is congruent to 2
modulo 3:

```jldoctest
julia> using DeRham, Oscar

julia> f7 = DeRham.fermat_hypersurface(3, 3, 7);

julia> FM7 = DeRham.frobenius_matrix_with_precision(f7, [1, 0]);

julia> H7 = DeRham.hasse_witt_matrix(FM7, 1);

julia> (rank(H7), DeRham.a_number(FM7, 1))
(1, 0)

julia> f5 = DeRham.fermat_hypersurface(3, 3, 5);

julia> FM5 = DeRham.frobenius_matrix_with_precision(f5, [1, 0]);

julia> H5 = DeRham.hasse_witt_matrix(FM5, 1);

julia> (rank(H5), DeRham.a_number(FM5, 1))
(0, 1)
```

## API

```@docs
DeRham.hasse_witt_matrix
DeRham.a_number
```
