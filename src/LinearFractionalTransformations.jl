# Linear fractional transformation "LFT"

module LinearFractionalTransformations
using StaticArrays

import Base: inv, show, hash, ==, *, isone
import Base: getindex


export LFT, call, RiemannInfinity, ∞


include("infinity.jl")


struct LFT{T}
    M::SMatrix{2,2,T,4}
    function LFT(a, b, c, d)
        if is_infinite_point(a) || is_infinite_point(b) ||
           is_infinite_point(c) || is_infinite_point(d)
            error("Arguments must be finite: " * string((a, b, c, d)))
        end

        if iszero(a * d - b * c)
            error("Singularity detected: " * string((a, b, c, d)))
        end
        a, b, c, d = promote(a, b, c, d)
        T = typeof(a)
        new{T}(SMatrix{2,2,T,4}(a, c, b, d))
    end
end

"""
    LFT
A linear fractional transformation `f(z)=(a*z+b)/(c*z+d)` where `a,b,c,d` are complex numbers.

Constructors:
* `LFT(a,b,c,d)`: create the function described above.
* `LFT()`: create the identity function `f(z) = z`. Equivalent to `LFT(1,0,0,1)`.
* `LFT(a,b,c)`: create the function with `a ↦ 0`, `b ↦ 1`, and `c ↦ ∞`.
* `LFT(a,aa,b,bb,c,cc)` or `LFT(a=>aa,b=>bb,c=>cc)`: creates the function with `a ↦ aa`, `b ↦ bb`, and `c ↦ cc`.
"""
function LFT(M::AbstractArray)
    return LFT(M[1, 1], M[1, 2], M[2, 1], M[2, 2])
end


function LFT(aaa::Pair, bbb::Pair, ccc::Pair)
    a, aa = aaa
    b, bb = bbb
    c, cc = ccc
    return LFT(a, aa, b, bb, c, cc)
end


function LFT()
    return LFT(1, 0, 0, 1)
end

# a --> 0, b-->1, and c--> oo
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

# a-->aa, b-->bb, c-->cc
function LFT(a, aa, b, bb, c, cc)
    f = LFT(a, b, c)
    g = LFT(aa, bb, cc)
    return inv(g) * f
end


#### Equality checking ####

# These need to be rewritten

"""
    isone(f::LFT)::Bool

Return `true` if `f` is the identity `LFT` and `false` otherwise.
"""
function isone(f::LFT)::Bool
    M = f.M
    return iszero(M[1, 2]) && iszero(M[2, 1]) && (M[1, 1] == M[2, 2])
end


function ==(f::LFT, g::LFT)::Bool
    return isone(f * inv(g))
end


# Inverse transformation
function inv(L::LFT{T})::LFT{T} where {T}
    a = L.M[1, 1]
    b = L.M[1, 2]
    c = L.M[2, 1]
    d = L.M[2, 2]
    return LFT(d, -b, -c, a)
end

# Composition
*(A::LFT, B::LFT) = LFT(A.M * B.M)

# Function application
function (A::LFT)(x)
    a = A.M[1, 1]
    b = A.M[1, 2]
    c = A.M[2, 1]
    d = A.M[2, 2]
    numer, denom = is_infinite_point(x) ? (a, c) : (a * x + b, c * x + d)
    return iszero(denom) ? infinity(denom) : numer * inv(denom)
end

# call(A::LFT, x::Number) = A[x]

function show(io::IO, L::LFT)
    print(
        io,
        "LFT( ",
        L.M[1, 1],
        " , ",
        L.M[1, 2],
        " , ",
        L.M[2, 1],
        " , ",
        L.M[2, 2],
        " )",
    )
end


function hash(f::LFT, h::UInt64 = UInt64(0))
    a = f.M[1, 1]
    b = f.M[1, 2]
    c = f.M[2, 1]
    d = f.M[2, 2]
    # The images at 0, 1 and ∞. `+ zero(x)` maps -0.0 to 0.0, which Julia hashes
    # apart, while a pole's infinity is returned untouched — it may be a non-numeric
    # ring (a Nemo field's), so it must not be added to anything.
    _image(x, y) = iszero(y) ? infinity(x) : (x * inv(y) + zero(x))
    return hash(_image(b, d), hash(_image(a + b, c + d), hash(_image(a, c), h)))
end

include("projection.jl")
include("deprecated.jl")



end # end of module "LFTs"
