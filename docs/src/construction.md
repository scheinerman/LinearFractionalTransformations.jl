```@meta
CurrentModule = LinearFractionalTransformations
```

# Construction

A linear fractional transformation is built from its four coefficients, from a `2×2`
matrix, or from the images of three points.

```julia
f = LFT(1, 2, 3, 4)             # z -> (z + 2) / (3z + 4)
g = LFT([1 2; 3 4])             # same map, from its matrix
h = LFT(0, 1, Inf, 1, 2, 3)     # 0 -> 1, 1 -> 2, Inf -> 3
```

The constructor overloads:

| Form | Result |
| --- | --- |
| `LFT(a, b, c, d)` | the map with those coefficients |
| `LFT(M::AbstractMatrix)` | from the `2×2` matrix `[a b; c d]` |
| `LFT(a, b, c)` | `a ↦ 0`, `b ↦ 1`, `c ↦ ∞` |
| `LFT(a, aa, b, bb, c, cc)` | `a ↦ aa`, `b ↦ bb`, `c ↦ cc` |
| `LFT(a => aa, b => bb, c => cc)` | pair form of the six-point constructor |
| `LFT()` | the identity, `LFT(1, 0, 0, 1)` |
| `LFT(::Type{T})` | the identity of type `T` |

The four coefficients are promoted to a common type, so `LFT(1, 2.0, 3, 4)` has
`Float64` coefficients. Infinite coefficients are rejected, as is a singular
coefficient vector with `ad - bc == 0`.

```@docs
LFT
```
