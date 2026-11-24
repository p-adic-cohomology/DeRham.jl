# Pairs of quadrics

Let $X=V(f,g)$ in $\mathbf P^{n-1}$ over a finite field of odd order $q$.
`quadric_pencil_point_count` computes its exact number of rational points,
including for singular intersections. Its cost includes $q+1$ ranks of
$n$-by-$n$ matrices. For odd $n$ and a smooth complete intersection,
`quadric_pencil_zeta` computes the entire rational zeta function by factoring
one polynomial of degree at most $n$. Both routines use Oscar finite-field
arithmetic. The full zeta routine uses no field enumeration or choice of
working precision.

```jldoctest
julia> using DeRham, Oscar

julia> _, x = polynomial_ring(GF(7), ["x$i" for i = 1:5]);

julia> f = sum(x .^ 2); g = sum((i-1)*x[i]^2 for i = 1:5);

julia> DeRham.quadric_pencil_point_count(f, g)
64

julia> Z = DeRham.quadric_pencil_zeta(f, g); P, T = polynomial_ring(ZZ, "T");

julia> Z == one(P) // ((1-T)*(1-7T)^4*(1+7T)^2*(1-49T))
true
```

## Counting identity

Write $M_{[s:t]}=sA+tB$ for the pencil of symmetric Gram matrices of $f,g$.
Above each ambient point, exactly one member vanishes unless the point lies on
$X$, in which case all $q+1$ vanish. Thus the sum of the projective quadric
counts equals $\#\mathbf P^{n-1}+q\#X$.

A member of rank $k$ contributes $\#\mathbf P^{n-2}$ plus a correction.
The correction vanishes for odd $k$. For even $k$ it equals
$\epsilon q^{n-1-k/2}$, where
$\epsilon=\chi((-1)^{k/2}\det M_{\rm nd})$ and $M_{\rm nd}$ is a
nonsingular quotient by the radical. Use the empty determinant $1$ for
rank zero. Indeed, diagonalize and apply additive-character orthogonality
to the affine zero count. Each nonzero diagonal coefficient $d$ contributes
$\chi(ad)G$ at frequency $a\ne 0$, with $G^2=\chi(-1)q$.
Summing over $a$, multiplying by the radical size, and projectivizing gives
the stated correction. The pivot columns of a row echelon form select a
coordinate complement to the kernel, so its principal minor computes this sign.
The incidence identity now gives $\#X=\#\mathbf P^{n-3}$ plus the sum of the
corrections divided by $q$.

## Recovering all extension counts

For a smooth pencil, the binary determinant has $n$ distinct geometric roots.
To see the criterion, work over the algebraic closure. A member of corank at
least two has a nonzero vector in its kernel on which a second pencil form
vanishes. At a corank-one root the determinant derivative is a nonzero scalar
times that second form on a generator of the kernel, by the adjugate formula.
These are exactly the common zeros with dependent gradients. The same argument
in the other affine chart covers infinity. An identically singular pencil
also supplies such a point by applying the same corank and derivative tests.
The code checks squarefreeness and allows a single simple root at infinity.

When $n$ is odd, the only corrections come from these roots, all of rank $n-1$.
A closed root of degree $e$ contributes over $\mathbf F_{q^a}$ precisely when
$e\mid a$. There it has $e$ conjugates and sign $\epsilon^{a/e}$, where
$\epsilon$ is the sign over its residue field. A nonzero signed principal
cofactor gives that sign. For a monic irreducible factor $h$, its resultant
with the cofactor is the norm to the base field, which preserves the quadratic
character. At infinity use a signed cofactor of $A$. Applying the exponential
definition of the zeta function to the counts over all extensions gives

```math
Z(X,T)^{-1}
 = \prod_{i=0}^{n-3}(1-q^i T)
   \prod_{\text{closed roots}}\bigl(1-\epsilon(q^{(n-3)/2}T)^e\bigr).
```

See Laura Mann Schueller,
[*The Zeta Function of a Pair of Quadratic Forms*](https://doi.org/10.4153/CMB-2001-024-7),
Canadian Mathematical Bulletin 44 (2001), 242–256.

## API

```@docs
DeRham.quadric_pencil_point_count
DeRham.quadric_pencil_zeta
```
