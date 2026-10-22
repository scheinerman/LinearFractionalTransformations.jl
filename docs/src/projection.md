```@meta
CurrentModule = LinearFractionalTransformations
```

# Projection

`stereo` maps between the complex plane and a unit sphere in three dimensions, and
`LFTQ` builds a transformation from a rotation of that sphere.

```julia
z = 1 + 2im
stereo(z)            # the point on the sphere
stereo(stereo(z))    # back to z, up to rounding

Q, _ = qr(randn(3, 3))     # a 3×3 orthogonal matrix
f = LFTQ(Matrix(Q))        # the LFT induced by Q
```

`LFTQ` requires a matrix `Q` with `Q*Q'` equal to the identity and `det(Q) == 1`.
`stereo` is its own inverse, so `stereo(stereo(z))` equals `z` up to floating-point
rounding.

```@docs
stereo
LFTQ
```
