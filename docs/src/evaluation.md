```@meta
CurrentModule = LinearFractionalTransformations
```

# Evaluation

An `LFT` is callable, mapping a point `z` to

```math
f(z) = \frac{a z + b}{c z + d}.
```

Evaluation at `Inf` is handled — as are `-Inf`, `complex(Inf)` and the sentinel
[`∞`](@ref) — and a pole (vanishing denominator) returns [`infinity`](@ref) of the
result field: that field's own infinity where it has one (`Inf` over `Float64`,
`complex(Inf)` over `ComplexF64`, `1//0` over ℚ, `Inf` over ℤ), and `∞` where it has
none (a Nemo number field). See [Infinity](infinity.md). Maps broadcast elementwise
over collections.

```julia
f = LFT(0.0, 1.0, 1.0, 0.0)      # z -> 1/z
f(0.0)                           # Inf   — a Float64, so the results stay concrete
f.([0.0, 2.0])                   # Vector{Float64}: works with sum, sort, …

LFT(0//1, 1//1, 1//1, 0//1)(0//1)   # 1//0 — ℚ's own infinity, still a Rational
```

```julia
f = LFT(1, 2, 3, 4)
f(0), f(1), f(Inf)
f.([0, 1, Inf])
```

```@docs
call
Base.broadcastable
```
