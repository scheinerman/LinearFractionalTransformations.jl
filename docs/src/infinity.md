```@meta
CurrentModule = LinearFractionalTransformations
```

# Infinity

A linear fractional transformation is a map on the Riemann sphere, so the point at
infinity plays two roles:

- **returned** — a map evaluated at a pole (vanishing denominator) returns it;
- **accepted** — `f(∞)`, `f(Inf)`, `f(-Inf)`, `f(complex(Inf))` are all recognised on
  input and evaluate to the limit `a/c` (and so is `1//0`, which is what a pole over
  ℚ returns — there, accepted and returned are the same value).

Which *value* stands for it depends on the field the division lands in — the field of
the denominator `c*z + d`. That is usually the coefficient field, but not always: over
`Int` coefficients a `Float64` input makes the whole computation `Float64`, and a pole
then returns the `Float64` infinity.

```julia
julia> fi = LFT(1, 2, 3, 5);   # integer coefficients

julia> fi(1)                   # …yet the computation is in Float64
0.375

julia> fi(-5/3)                # pole 3z + 5 = 0, reached through a Float64 input
Inf
```

The table is indexed by that field — the coefficient field in the common case where the
input matches it:

| field | a pole returns | why |
| --- | --- | --- |
| `Float64`, `Float32`, `BigFloat`, `Complex{Float64}`, … | `Inf`, `Inf32`, `big(Inf)`, `complex(Inf)`, … | the field has its own infinity, so the result keeps a concrete element type |
| `Rational` | `1//0` | ℚ's own division by zero: still a `Rational`, so the result stays in ℚ |
| `Integer` | `Inf` | ℤ has no infinity of its own, but ℤ's `/` already lands in `Float64` — a pole lands there too, matching the rest of the map |
| a Nemo number field, or any field type the package does not know | `∞` | no infinity could be named for it — only the field's own type knows what one looks like there |

`∞` appears only where the package cannot name the field's infinity. Extend `infinity`
and even that field supplies its own.

The point is *typing*: over `Float64` the promise that the result *is* a `Float64` keeps
`f.(grid)` a `Vector{Float64}` whose entries `sort`, `sum` and ordinary arithmetic work
on. The same holds over ℚ and ℤ — a pole does not knock the container down to `Any`:

```julia
julia> pole = LFT(0.0, 1.0, 1.0, 0.0);    # z -> 1/z

julia> pole(0.0)                          # the field's own infinity
Inf

julia> pole(0.0) == Inf
true

julia> exact = LFT(0//1, 1//1, 1//1, 0//1);   # the same map over ℚ

julia> exact(0//1)                        # ℚ's own infinity, still a Rational
1//0

julia> isinf(exact(0//1))                 # …and it answers to isinf
true

julia> eltype([exact(z) for z in (0//1, 1//1, 2//1)])
Rational{Int64}                           # the pole did not box the grid

julia> sort([exact(z) for z in (0//1, 1//1, 2//1)])
3-element Vector{Rational{Int64}}:
 1//2
 1//1
 1//0

julia> integer = LFT(1, 0, 2, 2);         # ℤ coefficients, pole at z = -1

julia> integer(0), integer(-1)            # what the rest of the map computes in…
(0.0, Inf)
```

`∞` hashes like `Inf` and compares equal to it (`∞ == Inf`), so the two are
interchangeable as `Dict`/`Set` keys. `1//0` does the same (`hash(1//0) == hash(Inf)`),
so a map hashes alike whichever of the three its own field supplies. `-Inf` stays
distinct, as it is in Julia.

The sentinel's arithmetic is deliberately minimal: `zero(∞)` and `∞ + ∞` are both `∞`, so
`∞` is its own additive identity and generic additive code (`sum`, `zeros`) works on it —
while `∞ + 1` is a `MethodError`, so arithmetic that runs into a pole fails loudly rather
than absorbing quietly. ℚ's own infinity cannot make that promise: `1//0 + 1` is `1//0`,
because Julia's rational arithmetic is total.

## Choosing the value a pole returns

[`infinity`](@ref) decides it, and it takes an *element* of the field, not a type: over
an exact field only an element carries its parent ring, which an extension needs to name
the field's own infinity.

```@docs
∞
RiemannInfinity
infinity
```

## Over a Nemo number field

```julia
using Nemo
R, x = polynomial_ring(QQ, "x")
K, a = number_field(x^2 + 1, "a")   # ℚ(i)

f = LFT(K(1), K(2), K(3), K(4))     # z -> (z + 2)/(3z + 4)

f(Nemo.inf)      # accepted on input (isinf is defined for it) → a/c = 1/3
f(-K(4)//K(3))   # the pole: 3·(-4/3) + 4 = 0 → ∞, since ℚ(i) has no infinity
```

`Nemo.inf` is AbstractAlgebra's positive infinity. Its type is not a `Real`, so it is
recognised through `isinf` instead. A pole over ℚ(i) returns `∞`, which is not a field
element, and `Nemo.inf` is accepted on the way in with no configuration.

To have a pole over ℚ(i) return something from the field, extend `infinity`:

```julia
LinearFractionalTransformations.infinity(x::AbsSimpleNumFieldElem) = infinity(parent(x))
```

## Known limits

**Pole detection is exact zero.** A pole is recognised by `iszero(denominator)`, so a
field whose zero is not decidable from a value never reaches that branch: over a ball
field (`AcbField`, `ArbField`), a denominator that *contains* zero is not `iszero`, and
the division proceeds, giving an indeterminate ball (`nan`) rather than `∞`. A different
pole test would be a different semantics, not a different representation.

**Hash agreement across coefficient types is best-effort.** The hash is computed from
the images `f(0)`, `f(1)`, `f(∞)` in their own field, and an exact value is not its
nearest float: `6/14 == 3//7` is `false`, so `LFT(1.0, 2.0, 3.0, 7.0)` and its ℚ twin
compare equal yet hash differently. Maps over different coefficient types are therefore
the same `Dict`/`Set` key only when their images agree exactly.

**A `Union` is not automatic for the fields that still use the sentinel.** Where a
field supplies its own infinity the container is that field, as above. Where none could
be named, the result is `Union{T, RiemannInfinity}`, and the widening hook decides
whether a container shows it: comprehensions and `map` go through `promote_typejoin`,
which this package extends to give that union, but a *literal* vector `[∞, x]` goes
through `promote_type`, which rejects a `Union`, so it stays a `Vector{Any}`.
Comprehensions and `map` are the case that matters for a grid; elsewhere, write
`Union{T, RiemannInfinity}[…]` to force it.
