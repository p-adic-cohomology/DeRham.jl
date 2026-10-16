# Artin–Mazur height

The Artin–Mazur height of a K3 surface is the height of its one-dimensional
formal Brauer group. It is an integer from 1 to 10, or infinity. Height 1 is
the ordinary case; infinite height is the supersingular case.

The polynomial form is for a smooth quartic K3 surface in `P^3`. The data form
takes the 22 coefficients of its primitive `H^2` L-polynomial over `F_q`;
it accepts raw `zeta_coefficients` output, `tate_twist` output, or
`normalized_tate_twist` output.

```@docs
DeRham.artin_mazur_height
```
