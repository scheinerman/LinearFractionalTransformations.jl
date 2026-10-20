# Infinity
#
# The coefficients may live in any field, including ones whose elements have no
# `isinf` and no floating-point infinity: an exact field such as ℚ, or a Nemo
# number field. A pole must then return an element of that field, or a sentinel
# when the field has no infinity to name.

"""
    RiemannInfinity

Singleton type whose unique instance `∞` is the point at infinity of the Riemann
sphere. It is the fallback: what a pole evaluates to when the field the division lands
in is one the package cannot name an infinity for — a Nemo number field, or a field
type it does not know. The fields it does know supply their own (`Inf` over `Float64`,
`1//0` over `Rational`), so ordinary work keeps a concrete element type.
"""
struct RiemannInfinity end

"""
    ∞

The point at infinity of the Riemann sphere: the one instance of
[`RiemannInfinity`](@ref), and the value a pole takes where the field has no infinity
of its own.
"""
const ∞ = RiemannInfinity()   # type \infty<TAB>

Base.show(io::IO, ::RiemannInfinity) = print(io, "∞")
# The `applicable` fallback below uses this to recognise `∞`.
Base.isinf(::RiemannInfinity) = true

# Just enough arithmetic for the sentinel to be its own additive identity, so that
# generic additive code (`sum`, `zeros`, `mapreduce` with an `init`) works.
# Numbers are deliberately left out: `∞ + 1` stays an error,
# so arithmetic that runs into a pole fails loudly instead of quietly absorbing.
Base.zero(::RiemannInfinity) = ∞
Base.:+(::RiemannInfinity, ::RiemannInfinity) = ∞

# One hash for every infinity: Julia already hashes `complex(Inf)` like `Inf`, so the
# whole pole family collapses onto one value.
Base.hash(::RiemannInfinity, h::UInt) = hash(Inf, h)

# `==` says what `hash` says, so the sentinel works as a key alongside the infinities
# it stands for (`Dict(Inf => 1)[∞]`). `-Inf` stays distinct, as it is in Julia.
# Restricted to `Number`, so an exotic infinity (AbstractAlgebra's `PosInf`) stays out
# until someone extends this for it.
Base.:(==)(::RiemannInfinity, x::Number) = (x == Inf)
Base.:(==)(x::Number, ::RiemannInfinity) = (x == Inf)

# A container that mixes `∞` with other values should be typed `Union` rather than
# `Any`: `collect` widens with `promote_typejoin`, whose default would give `Any`.
#
# The other side is left unconstrained because fields that fall back to `∞` — a Nemo
# number field, a user-defined one — are not `Number`s.
Base.promote_typejoin(::Type{T}, ::Type{RiemannInfinity}) where {T} = Union{T, RiemannInfinity}
Base.promote_typejoin(::Type{RiemannInfinity}, ::Type{T}) where {T} = Union{T, RiemannInfinity}
# Those two are ambiguous on their common point; this resolves it, and keeps `[∞, ∞]`
# a `Vector{RiemannInfinity}`.
Base.promote_typejoin(::Type{RiemannInfinity}, ::Type{RiemannInfinity}) = RiemannInfinity

# Fast path: `Real` and `Complex` define `isinf`, so this skips the guarded fallback
# below. Everything else takes the slow path.
is_infinite_point(z::Real) = isinf(z)
is_infinite_point(z::Complex) = isinf(z)
# The guard covers `∞`, AbstractAlgebra's `PosInf`/`NegInf` (neither `Real` nor
# `Complex`), and any type that defines `isinf`. Types with no `isinf` (Nemo
# number-field elements) are not infinities.
is_infinite_point(z) = applicable(isinf, z) && isinf(z)

"""
    infinity(x)

The value a pole evaluates to when the result of `m(z)` lives in the field of `x`:
that field's own infinity, wherever the package can name one — `Inf` over the
floating-point fields, `complex(Inf)` over `Complex`, `1//0` over `Rational`, and
`Inf` over `Integer` (whose `/` lands in `Float64` anyway).

[`RiemannInfinity`](@ref) (`∞`) is the fallback, for a field the package cannot name
an infinity for: a Nemo number field, or any field type it does not know.

It takes an *element*, not a type: over an exact field only an element carries its
parent ring, which is what an extension needs to name the field's own infinity. Extend
it for a new field type with `infinity(x::MyElem) = infinity(parent(x))`.
"""
infinity(x::T) where {T <: AbstractFloat} = oftype(x, Inf)
infinity(x::Integer) = one(x) / zero(x)    # ℤ has no infinity, but ℤ's `/` lands in Float64
infinity(x::Rational) = one(x) / zero(x)   # 1//0 — an infinity that stays a Rational
infinity(x::Complex) = complex(infinity(real(x)))
infinity(_) = ∞      # a field whose infinity could not be named: Nemo, or unknown

# The float infinity the package assumed before exact fields; kept for compatibility.
const complex_infinity = Inf + Inf * im
