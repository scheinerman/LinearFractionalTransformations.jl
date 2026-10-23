```@meta
CurrentModule = LinearFractionalTransformations
```

# LinearFractionalTransformations.jl

A linear fractional transformation is a rational function of the form

```math
f(z) = \frac{a z + b}{c z + d}, \qquad ad - bc \neq 0,
```

mapping the extended points of a field to itself. The package provides a lightweight
`LFT` type with constructors, evaluation, composition, and projective equality. The
coefficients may come from any field — `Int`, `Rational`, `Complex`, or an exact field
such as a Nemo number field.

## Installation

```julia
julia> ] add LinearFractionalTransformations
```

## Usage

```julia
using LinearFractionalTransformations

f = LFT(1, 2, 3, 4)               # z -> (z + 2) / (3z + 4)
f(0), f(1), f(Inf)                # evaluation, Inf included

g = LFT(0, 1, Inf, 1+im, 2+2im, 3+3im)   # 0 -> 1+i, 1 -> 2+2i, Inf -> 3+3i

h = LFT(0//1, 1//1, 1//1, 0//1)   # the same map z -> 1/z over ℚ
h(0//1)                           # over ℚ a pole returns 1//0, still a Rational
```

## Sections

- [Construction](construction.md) — build an `LFT`
- [Evaluation](evaluation.md) — apply a transformation to a point
- [Operations](operations.md) — invert, compose, compare, matrix, determinant, normalization
- [Infinity](infinity.md) — the point at infinity: what a pole returns, what is accepted
- [Projection](projection.md) — stereographic projection and transformations from rotations

## License

MIT.
