# LinearFractionalTransformations



This module defines a `LFT` data type to represent a complex *linear
fractional transformation*. This is a function on the extended
complex numbers (include complex infinity) defined by
```
f(z) = (a*z + b) / (c*z + d)
```
where `a,b,c,d` are (finite) complex numbers and `a*d-b*c != 0`.

These are also known as *Möbius transformations*.

## Constructors

The basic constructor takes four values:

```julia
julia> using LinearFractionalTransformations

julia> julia> f = LFT(1,2,3,4)
LFT( 1.0 + 0.0im , 2.0 + 0.0im , 3.0 + 0.0im , 4.0 + 0.0im )
```

Notice that the `LFT` is represented by a 2-by-2 complex matrix.
A `LFT` can also be defined by specifying a matrix.

```julia
julia> A = [1 2; 3 4];

julia> g = LFT(A)
LFT( 1.0 + 0.0im , 2.0 + 0.0im , 3.0 + 0.0im , 4.0 + 0.0im )
```

The identity `LFT` is constructed by `LFT()`:

```julia
julia> LFT()
LFT( 1.0 + 0.0im , 0.0 + 0.0im , 0.0 + 0.0im , 1.0 + 0.0im )
```

Given (complex) numbers `a,b,c` (including `Inf`) we can construct
a `LFT` that maps `a` to 0, `b` to 1, and `c` to infinity.

```julia
julia> f = LFT(2,5,-1)
LFT( 6.0 + 0.0im , -12.0 + 0.0im , 3.0 + 0.0im , 3.0 + 0.0im )

julia> f[2]
0.0 + 0.0im

julia> f[5]
1.0 + 0.0im

julia> f[-1]
Inf + Inf*im
```

Finally, we provide a constructor for mapping a given triple of values
`(a,b,c)` to another triple `(aa,bb,cc)`. The syntax is
`LFT(a,aa,b,bb,c,cc)` or `LFT(a=>aa, b=>bb, c=>cc)`. Here's an example:

```julia
julia> f = LFT(1,2+im, 3,Inf, 4,1-im)
LFT( 5.0 + 1.0im , -17.0 - 7.0im , 3.0 + 0.0im , -9.0 + 0.0im )

julia> f[1]
2.0 + 1.0im

julia> f[3]
Inf + Inf*im

julia> f[4]
1.0 - 1.0im
```


#### Under the hood

The matrix representing a `LFT` object is held in a field named `:M`.

```julia
julia> f = LFT(1,2,3)
LFT( -1.0 + 0.0im , 1.0 + 0.0im , 1.0 + 0.0im , -3.0 + 0.0im )

julia> f.M
2x2 Array{Complex{Float64},2}:
 -1.0+0.0im   1.0+0.0im
  1.0+0.0im  -3.0+0.0im
```

## Operations

### Function application

Since a `LFT` is a function, the most basic operation we may wish to
perform is applying that function of a complex number. That's done
with  `f(x)`:

```julia
julia> f = LFT(3,2,1,1)
LFT( 3.0 + 0.0im , 2.0 + 0.0im , 1.0 + 0.0im , 1.0 + 0.0im )

julia> f(1)
2.5 + 0.0im

julia> f(0)
2.0 + 0.0im

julia> f(-1)
Inf + Inf*im

julia> f(Inf)
3.0 + 0.0im

julia> f(1+2im)
2.75 + 0.25im
```


### Composition and inverse

The `*` operation is used for function composition.

```julia
julia> f = LFT(3,2,1,1);

julia> g = LFT(0,1,-1,2);

julia> f*g
LFT( -2.0 + 0.0im , 7.0 + 0.0im , -1.0 + 0.0im , 3.0 + 0.0im )

julia> g*f
LFT( 1.0 + 0.0im , 1.0 + 0.0im , -1.0 + 0.0im , 0.0 + 0.0im )
```


The inverse of a `LFT` is computed with `inv`:

```julia
julia> f = LFT(1,2,3,4);

julia> g = inv(f)
LFT( 4.0 + 0.0im , -2.0 - 0.0im , -3.0 - 0.0im , 1.0 + 0.0im )

julia> f*g
LFT( -2.0 + 0.0im , 0.0 + 0.0im , 0.0 + 0.0im , -2.0 + 0.0im )
```

Notice that the matrix representing `f*g` is a scaled version of the
identity matrix.

## Equality checking

We can use `==` or `isequal` to check if two `LFT` objects are
equal. Note that there is no unique matrix representation for a `LFT`
object and we might have that `f` and `g` are equal, but `f.M` and
`g.M` are different.

```julia
julia> f = LFT(1,2,3,4);

julia> g = LFT(-2,-4,-6,-8);

julia> f==g
true

julia> f.M == g.M
false
```


## Stereographic Projection and Linear Fractional Transformations

### The `stereo` function

The function `stereo` maps points in the complex plane to points on a unit three-dimensional sphere centered at `[0,0,0]` via [sterographic projection](https://en.wikipedia.org/wiki/Stereographic_projection).

The north pole, `[0,0,1]`, corresponds to complex infinity and the 
south pole, `[0,0,-1]`, corresponds to `0+0im`.

For a complex number `z`, `stereo(z)` maps `z` to the sphere. This may also be invoked as 
`stereo(x,y)`.

For a (unit) three-dimensional vector `v`, `stereo(v)` returns the complex number by projecting 
`v` to the complex plane. This may also be invoked as `stereo(x,y,z)`. Note that the function does not check if `v` is a unit vector.

Note that `stereo` is its own inverse. 
That is, for a complex number `v`, we have `stereo(stereo(v))` should equal `v` were it not for 
possible roundoff errors. Likewise for a unit three-dimensional real vector.
```julia
julia> z = 3-4im
3 - 4im

julia> stereo(stereo(z))
3.000000000000002 - 4.000000000000003im

julia> u = [1,2,2]/3       # this is a unit vector
3-element Vector{Float64}:
 0.3333333333333333
 0.6666666666666666
 0.6666666666666666

julia> stereo(stereo(u))
3-element Vector{Float64}:
 0.3333333333333333
 0.6666666666666666
 0.6666666666666666
```




### Creating a `LFT` from a unitary matrix: `LFTQ`


Linear fractional transformations may be considered a mapping of a complex point to the unit sphere, followed by a rotation of the sphere, followed by a projection back to the complex plane. 

Specifically, if `Q` is a real, orthogonal 3-by-3 matrix [so `Q*Q'` is the identity and `det(Q)` equals 1, i.e., `Q` is in [SO(3)](https://en.wikipedia.org/wiki/Stereographic_projection)], the function `LFTQ(Q)` returns a linear fractional transformation `F` with the property that for complex `z`, we have `F(z)` equal to `stereo(Q*(stereo(v)))`. Roundoff errors may occur.
```julia
julia> using LinearAlgebra, LinearFractionalTransformations

julia> Q,R = qr(randn(3,3));  # create a randon orthogonal matrix

julia> det(Q)
1.0

julia> Q' == inv(Q)   # verify that Q is in SO(3)
true

julia> v = 5-im
5 - 1im

julia> F = LFTQ(Q);

julia> F(v)
1.0556981607448988 + 1.5029664004078547im

julia> stereo(Q*stereo(v))
1.055698160744897 + 1.5029664004078531im   #  slightly different from F(v)
```

