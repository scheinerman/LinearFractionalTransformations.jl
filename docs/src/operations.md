```@meta
CurrentModule = LinearFractionalTransformations
```

# Operations

Transformations compose, invert, and compare as projective maps — equal up to a
scalar multiple of the coefficients.

```julia
g = LFT(0, 1, 1, 0)                        # z -> 1/z
isone(g * inv(g))                          # true
LFT(1, 2, 3, 4) == LFT(2, 4, 6, 8)         # scalar multiples are equal
```

Composition makes the maps a monoid, so they have powers and an identity too — the
identity of `f`'s *own* coefficient type, not of a default one.

```julia
f = LFT(1, 2, 3, 4)
f^3 == f * f * f                   # powers compose, by squaring
f^0 == one(f)                      # the identity of f's coefficient type
f^-1 == inv(f)                     # a negative power is a power of the inverse
```

```@docs
Base.inv
Base.:(*)
Base.:∘
Base.one
Base.:(^)
Base.isone
Base.:(==)
Base.eltype
Base.copy
Base.hash
```

## Linear algebra

A transformation carries a `2×2` coefficient matrix `[a b; c d]`, with the
usual determinant and projective normalization.

```julia
f = LFT(1, 2, 3, 4)
Matrix(f)          # [1 2; 3 4]
det(f)             # 1*4 - 2*3
normalize(f)       # det == 1 (needs sqrt in the coefficient field)
```

```@docs
Base.Matrix
det
normalize
```
