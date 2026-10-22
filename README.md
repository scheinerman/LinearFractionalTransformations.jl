# LinearFractionalTransformations

## New in Version 0.3

The coefficients may come from any field, not just `Complex{Float64}`. 
See *Other coefficient fields* below. 

## Overview


This module defines an `LFT` data type to represent a *linear fractional transformation*,
a function on the extended points of a field defined by
```
f(z) = (a*z + b) / (c*z + d)
```
where `a,b,c,d` are finite elements of a field and `a*d-b*c != 0`.

These are also known as *Möbius transformations*.

## Constructors

The basic constructor takes four coefficients:

```
julia> using LinearFractionalTransformations

julia> f = LFT(1,2,3,4)
LFT( 1 , 2 , 3 , 4 )
```

The coefficients are stored in four fields and promoted to a common type. `Matrix(f)`
gives the 2-by-2 coefficient matrix:

```
julia> Matrix(f)
2×2 Matrix{Int64}:
 1  2
 3  4
```

An `LFT` can also be defined by specifying a matrix:

```
julia> A = [1 2; 3 4];

julia> g = LFT(A)
LFT( 1 , 2 , 3 , 4 )
```

The identity `LFT` is constructed by `LFT()`, or for a given coefficient type by
`LFT(T)`:

```
julia> LFT()
LFT( 1 , 0 , 0 , 1 )

julia> LFT(ComplexF64)
LFT( 1.0+0.0im , 0.0+0.0im , 0.0+0.0im , 1.0+0.0im )
```

Given (complex) numbers `a,b,c` (including `Inf` or `∞`, see
[Other coefficient fields](#other-coefficient-fields)) we can construct
a `LFT` that maps `a` to 0, `b` to 1, and `c` to infinity.

```
julia> f = LFT(2,5,-1)
LFT( 6 , -12 , 3 , 3 )

julia> f[2]
0.0

julia> f[5]
1.0

julia> f[-1]
Inf
```

Finally, we provide a constructor for mapping a given triple of values
`(a,b,c)` to another triple `(aa,bb,cc)`. The syntax is
`LFT(a,aa,b,bb,c,cc)` or `LFT(a=>aa, b=>bb, c=>cc)`. Here's an example:

```
julia> f = LFT(1,2+im, 3,Inf, 4,1-im)
LFT( 5+1im , -17-7im , 3+0im , -9+0im )

julia> f[1]
2.0 + 1.0im

julia> f[3]
Inf + 0.0im

julia> f[4]
1.0 - 1.0im
```


#### Under the hood

The coefficients are four fields of the `LFT` object:

```
julia> f = LFT(1,2,3)

julia> f.a, f.b, f.c, f.d
(-1, 1, 1, -3)

julia> Matrix(f)
2×2 Matrix{Int64}:
 -1   1
  1  -3
```

`f.M` returns that same matrix, for code written against earlier versions.

## Operations

### Function application

Since a `LFT` is a function, the most basic operation we may wish to
perform is applying that function to a point. That's done
with  `f(x)`:

```
julia> f = LFT(3,2,1,1)
LFT( 3 , 2 , 1 , 1 )

julia> f(1)
2.5

julia> f(0)
2.0

julia> f(-1)
Inf

julia> f(Inf)
3.0

julia> f(1+2im)
2.75 + 0.25im
```

Over a field without a nameable infinity, a pole returns `∞` instead of `Inf`; see
[Other coefficient fields](#other-coefficient-fields).


### Composition and inverse

The `*` operation is used for function composition, and `∘` is its synonym.

```
julia> f = LFT(3,2,1,1);

julia> g = LFT(0,1,-1,2);

julia> f*g
LFT( -2 , 7 , -1 , 3 )

julia> g*f
LFT( 1 , 1 , -1 , 0 )
```


The inverse of a `LFT` is computed with `inv`:

```
julia> f = LFT(1,2,3,4);

julia> g = inv(f)
LFT( 4 , -2 , -3 , 1 )

julia> f*g
LFT( -2 , 0 , 0 , -2 )
```

Notice that the matrix representing `f*g` is a scaled version of the
identity matrix.

### Identity, determinant and normalization

`one(f)` is the identity of `f`'s coefficient type, `det(f)` is the determinant of the
coefficient matrix, and `normalize(f)` is the same map with determinant 1.

```
julia> f = LFT(1,2,3);

julia> one(f)
LFT( 1 , 0 , 0 , 1 )

julia> isone(one(f))
true

julia> det(f)
2

julia> normalize(LFT(4.0,0,0,1))
LFT( 2.0 , 0.0 , 0.0 , 0.5 )

julia> det(normalize(LFT(4.0,0,0,1)))
1.0
```

`normalize` requires `sqrt` in the coefficient field.

### Powers

`f^k` is the `k`-fold composite of `f` with itself, computed by repeated squaring. A
negative exponent uses the inverse, and `f^0` is `one(f)`.

```
julia> f = LFT(1,2,3);

julia> f^3
LFT( -6 , 14 , 14 , -34 )

julia> f^3 == f*f*f
true
```

### Scalar multiples

Multiplying by a scalar scales every coefficient. The map is unchanged.

```
julia> f = LFT(1,2,3);

julia> 2 * f
LFT( -2 , 2 , 2 , -6 )

julia> (2 * f) == f
true
```

### Broadcasting

`f` broadcasts over a collection of points, so `f.(points)` applies it elementwise.

```
julia> f = LFT(2,0,0,1);    # f(z) = 2z

julia> f.([0, 1, 2])
3-element Vector{Float64}:
 0.0
 2.0
 4.0
```

`eltype(f)` is the coefficient type, and `copy(f)` returns `f` unchanged.

## Equality checking

We can use `==` or `isequal` to check if two `LFT` objects are
equal. Note that there is no unique matrix representation for a `LFT`
object and we might have that `f` and `g` are equal, but `Matrix(f)` and
`Matrix(g)` are different.

```
julia> f = LFT(1,2,3,4);

julia> g = LFT(-2,-4,-6,-8);

julia> f==g
true

julia> Matrix(f) == Matrix(g)
false
```


## Other coefficient fields

The coefficients may come from any field, not only `Complex{Float64}`. The
constructors therefore do not require `Number` arguments, and a pole — a point
`z` where `c*z + d == 0` — evaluates to the coefficient field's own infinity.

Over ℚ, for example, a pole is `1//0`, and the value of `f(z)` stays exact:

```
julia> f = LFT(0//1, 1//1, 1//1, 0//1)    # z ↦ 1/z over ℚ
LFT( 0 , 1 , 1 , 0 )

julia> f(2//1) == 1//2
true

julia> f(0//1) == 1//0
true
```

Over the floating-point types a pole is `Inf` (or `Inf32` for `Float32`), and
over `Complex` it is `Inf + 0im`.

Some fields have no infinity to name: an element of a number field, such as those
of Nemo, has neither `isinf` nor an infinite value. A pole then returns the
exported sentinel `∞`, an instance of `RiemannInfinity`:

```
julia> ∞
∞

julia> isinf(∞)
true
```

`∞` may be used wherever `Inf` is. The three-point constructor accepts it as one
of the points:

```
julia> f = LFT(∞, 0.0, 1.0);    # ∞ ↦ 0, 0 ↦ 1, 1 ↦ ∞

julia> f(∞)
0.0

julia> f(1.0)
Inf
```

It is not a coefficient, though — the four-argument constructor rejects it, as it
rejects `Inf`:

```
julia> LFT(∞, 0, 0, 1)
ERROR: Arguments must be finite: (∞, 0, 0, 1)
```

If you work with a field type of your own that can name an infinity, extend
`LinearFractionalTransformations.infinity` for one of its elements, and a pole
returns that value instead of `∞`.

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
```
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
```
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

