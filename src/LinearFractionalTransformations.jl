# Linear fractional transformation "LFT"

module LinearFractionalTransformations

import Base: inv, show, hash, ==, *, isone, getindex

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
    LFT()

The identity, `LFT(1, 0, 0, 1)`.
"""
LFT() = LFT(1, 0, 0, 1)

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
    if a == b || b == c || a == c
        error("Three arguments must be distinct: " * string((a, b, c)))
    end

    if is_infinite_point(a)
        return LFT(0, b - c, 1, -c)
    end

    if is_infinite_point(b)
        return LFT(1, -a, 1, -c)
    end

    if is_infinite_point(c)
        return LFT(1, -a, 0, b - a)
    end

    # if all args are finite
    aa = (b - c)
    bb = (-a) * (b - c)
    cc = (b - a)
    dd = (-c) * (b - a)
    return LFT(aa, bb, cc, dd)
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
    ==(m::LFT, n::LFT)

Return `true` if `m` and `n` are the same projective transformation, i.e. if
`m * inv(n)` is the identity.
"""
==(m::LFT, n::LFT)::Bool = isone(m * inv(n))

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
    Matrix(m::LFT)

The coefficient matrix `M = [a b; c d]`.
"""
function Base.Matrix(m::LFT)
    (; a, b, c, d) = m
    return [a b; c d]
end


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
