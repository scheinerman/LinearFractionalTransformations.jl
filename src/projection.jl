# Functions related to sterographic projection to/from a three-dimensional
# unit sphere centered at (0,0) = 0+0im

"""
    stereo 
Stereographic projection between the complex plane and a three-dimensional unit sphere 
centered at `[0,0,0]`. The north pole, `[0,0,1]`, corresponds to complex infinity and the 
south pole, `[0,0,-1]`, corresponds to `0+0im`.

For a complex number `z`, `stereo(z)` maps `z` to the sphere. This may also be invoked as 
`stereo(x,y)`.

For a (unit) three-dimensional vector `v`, `stereo(v)` returns the complex number by projecting 
`v` to the complex plane. This may also be invoked as `stereo(x,y,z)`.
"""
function stereo(z::Number)::Vector{Real}
    if isinf(z)
        return [0, 0, 1]
    end

    X, Y = reim(z)

    x = 2X
    y = 2Y
    z = X^2 + Y^2 - 1

    d = 1 / (1 + X^2 + Y^2)

    return d * [x, y, z]
end



function stereo(X::Real, Y::Real)::Vector
    if isinf(X) || isinf(Y)
        return [0, 0, 1]
    end

    x = 2X
    y = 2Y
    z = X^2 + Y^2 - 1

    d = 1 / (1 + X^2 + Y^2)

    return d * [x, y, z]
end


stereo(z::Complex) = stereo(reim(z)...)


function stereo(x::Real, y::Real, z::Real)::Complex
    if z == 1
        return Inf + im * Inf
    end
    X = x / (1 - z)
    Y = y / (1 - z)
    return X + im * Y
end


function stereo(v::Vector{T})::Complex where {T<:Real}
    stereo(v...)
end


"""
    LFTQ
Create a linear fractional transformation from a 3-by-3 unitary matrix. 

Given a 3-by-3 real matrix `Q` with `Q*Q'` equal to the identity and `det(Q)` equal to `1`,
create a `LFT` that maps a complex number `v` to
`stereo(Q*stereo(v))`.
"""
function LFTQ(Q::AbstractMatrix)::LFT

    zz = [0 + 0im, 1 + 0im, 0 + im]

    uu = stereo.(zz)
    uu = [Q * u for u in uu]

    ww = [stereo(u) for u in uu]

    LFT(zz[1], ww[1], zz[2], ww[2], zz[3], ww[3])

end

export stereo, LFTQ