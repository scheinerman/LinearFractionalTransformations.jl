# Linear fractional transformation "LFT"

module LinearFractionalTransformations

import Base: inv, show, hash, ==, *, isone, getindex
import LinearAlgebra: det, normalize

include("infinity.jl")

export LFT, call, RiemannInfinity, ∞


"""
    LFT{T} <: Function

A linear fractional transformation `f(z) = (a*z + b) / (c*z + d)` with coefficients
`a, b, c, d` in the field `T` and `a*d - b*c != 0`.

The coefficients are four fields of the same name, promoted to a common type. That
type may be any field: `Int`, `Rational`, `Complex`, or an exact field such as a Nemo
number field.

Constructors:
* `LFT(a, b, c, d)`: the map with those coefficients. Infinite and singular
  coefficient vectors are rejected.
* `LFT(M)`: the map with coefficient matrix `M = [a b; c d]`.
* `LFT()`: the identity, `LFT(1, 0, 0, 1)`.
* `LFT(a, b, c)`: the map `a ↦ 0`, `b ↦ 1`, `c ↦ ∞`. `Inf` and `∞` are permitted.
* `LFT(a, aa, b, bb, c, cc)` or `LFT(a=>aa, b=>bb, c=>cc)`: the map `a ↦ aa`,
  `b ↦ bb`, `c ↦ cc`.
"""
struct LFT{T} <: Function
    a::T
    b::T
    c::T
    d::T

    function LFT(a, b, c, d)
        if is_infinite_point(a) || is_infinite_point(b) ||
           is_infinite_point(c) || is_infinite_point(d)
            error("Arguments must be finite: " * string((a, b, c, d)))
        end
        if iszero(a * d - b * c)
            error("Singularity detected: " * string((a, b, c, d)))
        end
        a, b, c, d = promote(a, b, c, d)
        new{typeof(a)}(a, b, c, d)
    end
end
# `f.M` was the `SMatrix` of coefficients. It is still the same matrix, built on demand.
Base.getproperty(m::LFT, s::Symbol) = s === :M ? Matrix(m) : getfield(m, s)
Base.propertynames(::LFT) = (:a, :b, :c, :d, :M)


"""
    LFT(M::AbstractArray)

The transformation whose coefficient matrix is `M = [a b; c d]`.
"""
LFT(M::AbstractArray) = LFT(M[1, 1], M[1, 2], M[2, 1], M[2, 2])

"""
    LFT(::Type{T}=Int)

Return the identity LFT transformation of type `T`.
"""
LFT(::Type{T}=Int) where T = LFT(one(T), zero(T), zero(T), one(T))

"""
    LFT(a => aa, b => bb, c => cc)

The transformation sending `a ↦ aa`, `b ↦ bb`, `c ↦ cc`.
"""
function LFT(aaa::Pair, bbb::Pair, ccc::Pair)
    a, aa = aaa
    b, bb = bbb
    c, cc = ccc
    return LFT(a, aa, b, bb, c, cc)
end

"""
    LFT(a, b, c)

The transformation sending `a ↦ 0`, `b ↦ 1`, `c ↦ ∞`. Values of `Inf` or `∞` are
permitted.
"""
function LFT(a, b, c)
    if isequal(a, b) || isequal(b, c) || isequal(a, c)
        throw(ArgumentError("Three arguments must be distinct: $((a, b, c))"))
    end

    # Build 0 and 1 from a finite argument, so an exact field keeps its own type.
    if is_infinite_point(a)
        return LFT(zero(b), b - c, one(b), -c)
    elseif is_infinite_point(b)
        return LFT(one(a), -a, one(a), -c)
    elseif is_infinite_point(c)
        return LFT(one(a), -a, zero(a), b - a)
    end

    return LFT(b - c, -a * (b - c), b - a, -c * (b - a))
end

"""
    LFT(a, aa, b, bb, c, cc)

The transformation sending `a ↦ aa`, `b ↦ bb`, `c ↦ cc`. Values of `Inf` or `∞` are
permitted.
"""
LFT(a, aa, b, bb, c, cc) = inv(LFT(aa, bb, cc)) * LFT(a, b, c)
# `a ↦ 0`, `b ↦ 1`, `c ↦ ∞`
# `0 ↦ aa`, `1 ↦ bb`, `∞ ↦ cc`


#### Equality and identity ####

"""
    isone(m::LFT)

Return `true` if `m` is the identity transformation.
"""
function isone(m::LFT)::Bool
    (; a, b, c, d) = m
    return iszero(b) && iszero(c) && iszero(a-d)
end

"""
    one(m::LFT)

The identity of `m`'s coefficient type, built from its coefficients, so that
`one(m) * m == m`.
"""
function Base.one(m::LFT)
    (; a, b, c, d) = m
    return LFT(one(a), zero(b), zero(c), one(d))
end

"""
    ==(m::LFT, n::LFT)

Return `true` if `m` and `n` are the same projective transformation, i.e. if
`m * inv(n)` is the identity.
"""
==(m::LFT, n::LFT)::Bool = isone(m * inv(n))

Base.eltype(::LFT{T}) where T = T

"""
    copy(m::LFT)

Return `m`: the coefficients are immutable, so there is nothing to copy.
"""
Base.copy(m::LFT) = m

"""
    hash(m::LFT, h::UInt64 = UInt64(0))

Hash from the images of `0`, `1` and `∞` — `b/d`, `(a+b)/(c+d)` and `a/c` — which
determine `m` and are invariant under scaling its coefficients. A pole is hashed as
the field's own infinity, so scalar multiples hash alike over any field.
"""
function hash(m::LFT, h::UInt64 = UInt64(0))
    (; a, b, c, d) = m
    _image(x, y) = iszero(y) ? infinity(x) : (x * inv(y) + zero(x))
    Z = _image(b, d)       # m(0)
    U = _image(a+b, c+d)   # m(1)
    I = _image(a, c)       # m(∞)
    return hash(Z, hash(U, hash(I, h)))
end


#### Operations ####

"""
    inv(m::LFT)

The inverse transformation, with coefficient matrix `[d -b; -c a]`.
"""
function inv(m::LFT)
    (; a, b, c, d) = m
    return LFT(d, -b, -c, a)
end

"""
    *(m::LFT, n::LFT)

Compose the transformations: `(m * n)(z) == m(n(z))`.
"""
function *(m::LFT, n::LFT)
    (; a, b, c, d) = n
    e, f, g, h = a, b, c, d
    (; a, b, c, d) = m
    return LFT(a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)
end

"""
    *(λ, m::LFT)
    *(m::LFT, λ)

Scale every coefficient of `m` by `λ`. The map is unchanged.
"""
function *(λ, m::LFT)
    (; a, b, c, d) = m
    return LFT(λ * a, λ * b, λ * c, λ * d)
end
*(m::LFT, λ) = *(λ, m)

Base.:∘(m::LFT, n::LFT) = m * n

"""
    m^k::LFT

The `k`-fold composite of `m` with itself, by repeated squaring: `m^0` is `one(m)`,
and a negative exponent falls to `inv(m)^(-k)`.
"""
function Base.:^(m::LFT, k::Integer)
    k < 0 && return inv(m)^(-k)
    p, square = one(m), m
    while k > 0
        isodd(k) && (p = p * square)
        square = square * square
        k >>= 1
    end
    return p
end

"""
    Matrix(m::LFT)

The coefficient matrix `M = [a b; c d]`.
"""
function Base.Matrix(m::LFT)
    (; a, b, c, d) = m
    return [a b; c d]
end

"""
    det(m::LFT)

The determinant `a*d - b*c` of the coefficient matrix.
"""
function det(m::LFT)
    (; a, b, c, d) = m
    return a * d - b * c
end

"""
    normalize(m::LFT)

Return a transformation equal to `m` whose coefficient matrix has determinant 1.
Requires `sqrt` in the coefficient field.
"""
normalize(m::LFT) = m * inv(sqrt(det(m)))

"""
    broadcastable(m::LFT)

Return `Ref(m)`, so that `m.(z)` broadcasts elementwise rather than treating `m` as a
collection.
"""
Base.broadcastable(m::LFT) = Ref(m)


#### Evaluation ####

"""
    (m::LFT)(z)

Apply `m`: `m(z) = (a*z + b) / (c*z + d)`.

`Inf`, `∞` and any value whose type defines `isinf` evaluate to the limit `a/c`. A
pole evaluates to the coefficient field's own infinity — `Inf` over `Float64`,
`complex(Inf)` over `Complex`, `1//0` over `Rational` — or to `∞` where the field has
none to name.
"""
function (m::LFT)(z)
    (; a, b, c, d) = m
    numer, denom = is_infinite_point(z) ? (a, c) : (a * z + b, c * z + d)
    return iszero(denom) ? infinity(denom) : numer * inv(denom)
end

"""
    call(m::LFT, z)

Apply `m` to `z`; the same as `m(z)`.
"""
call(m::LFT, z) = m(z)


#### Display ####

function show(io::IO, m::LFT)
    (; a, b, c, d) = m
    print(IOContext(io, :compact => true), "LFT( ", a, " , ", b, " , ", c, " , ", d, " )")
end

# Without this, the REPL prints Base's generic-function form: `LFT <: Function`.
show(io::IO, ::MIME"text/plain", m::LFT) = show(io, m)


include("projection.jl")
include("deprecated.jl")

end # end of module "LFTs"
