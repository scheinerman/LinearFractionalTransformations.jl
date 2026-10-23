using Test, LinearFractionalTransformations, LinearAlgebra


@testset "Basics" begin
    f = LFT(1, 0, 0, 1)  # identity
    @test f(1 + im) == 1 + im


    f = LFT(2, 0, 0, 1)  # f(x) = 2x
    @test f(-im) == -2im

    g = inv(f)
    @test g(im) == 0.5im

    g = f * f
    @test g(2 - im) == 4 * (2 - im)

    f = LFT(1, 2 + im, 3, 0, 1 - im, 5)
    @test f(1) == 2 + im
    @test f(3) == 0
    @test f(1 - im) == 5

    g = inv(f)
    @test g(0) == 3

    @test LFT() == g * f

    f = LFT(0, 1, 1, 0)
    @test f(Inf) == 0
    @test isinf(f(0))


    f = LFT(3 - im, 2 + im, 6)
    @test f(3 - im) == 0
    @test f(2 + im) == 1
    @test isinf(f(6))

    f = LFT(3 => im, 0 => -5, im => -2)
    @test f(3) == im
    @test f(0) == -5
    @test f(im) == -2

end


# check out LFTQ 
@testset "LFTQ" begin
    for k = 1:5
        Q = [1 0; 0 -1]
        while det(Q) < 0.5
            Q, R = qr(randn(3, 3))
        end
        F = LFTQ(Matrix(Q))
        v = randn() + randn() * im

        w1 = stereo(Q * stereo(v))
        w2 = F(v)
        @test abs(w1 - w2) <= 1e-8

        z = stereo(Q * stereo(0))
        @test abs(z - F(0)) <= 1e-8

        z = stereo(Q * stereo(1))
        @test abs(z - F(1)) <= 1e-8

        z = stereo(Q * stereo(Inf))
        @test abs(z - F(Inf)) <= 1e-8

    end
end


@testset "LFT as a function" begin
    f = LFT(1, 2, 3, 4)
    @test f isa Function
    @test (f.a, f.b, f.c, f.d) == (1, 2, 3, 4)
    @test eltype(f) == Int
    @test copy(f) === f
    @test call(f, 1) == f(1)

    @test Matrix(f) == [1 2; 3 4]
    @test f.M == [1 2; 3 4]
    @test repr("text/plain", f) == "LFT( 1 , 2 , 3 , 4 )"   # REPL display, not the Function form
    @test LFT(f.M) == f
    @test det(f) == 1 * 4 - 2 * 3

    g = normalize(LFT(4.0, 0.0, 0.0, 1.0))
    @test g == LFT(4.0, 0.0, 0.0, 1.0)
    @test det(g) == 1.0
    @test g.a == 2.0 && g.d == 0.5

    @test isone(one(f))
    @test one(f) * f == f
    @test f^0 == one(f)
    @test f^5 == reduce(*, fill(f, 5))
    @test f^-1 == inv(f)
    @test inv(f) * f == one(f)

    @test (f ∘ f) == f * f
    @test 2 * f == f * 2
    @test (2 * f)(1 + im) == f(1 + im)

    @test f.([0, 1]) == [f(0), f(1)]
end


# A field element with no `isinf` — the shape of a Nemo number-field element.
struct ToyExact
    x::Rational{Int}
end
Base.:+(u::ToyExact, v::ToyExact) = ToyExact(u.x + v.x)
Base.:-(u::ToyExact, v::ToyExact) = ToyExact(u.x - v.x)
Base.:-(u::ToyExact) = ToyExact(-u.x)
Base.:*(u::ToyExact, v::ToyExact) = ToyExact(u.x * v.x)
Base.inv(u::ToyExact) = ToyExact(inv(u.x))
Base.:(==)(u::ToyExact, v::ToyExact) = u.x == v.x
Base.hash(u::ToyExact, h::UInt) = hash(u.x, h)
Base.iszero(u::ToyExact) = iszero(u.x)
Base.zero(::Type{ToyExact}) = ToyExact(0 // 1)
Base.one(::Type{ToyExact}) = ToyExact(1 // 1)
Base.zero(u::ToyExact) = ToyExact(0 // 1)


@testset "exact fields and infinity" begin
    @testset "the field's own infinity" begin
        f = LFT(0 // 1, 1 // 1, 1 // 1, 0 // 1)          # z -> 1/z over ℚ
        @test f(0 // 1) === 1 // 0
        @test f(1 // 1) === 1 // 1
        @test f(2 // 1) === 1 // 2
        @test map(f, [0 // 1, 2 // 1]) isa Vector{Rational{Int}}

        c = Complex{Rational{Int}}
        g = LFT(c(0), c(1), c(1), c(0))
        @test g(c(0)) === complex(1 // 0)
        @test g(c(1)) === c(1)
        @test map(g, [c(0), c(1)]) isa Vector{Complex{Rational{Int}}}
        @test [g(t) for t in (c(0), c(1))] isa Vector{Complex{Rational{Int}}}

        @test LFT(0.0, 1.0, 1.0, 0.0)(0.0) === Inf
        @test LFT(0.0f0, 1.0f0, 1.0f0, 0.0f0)(0.0f0) === Inf32
        @test LFT(0.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 0.0 + 0.0im)(0.0 + 0.0im) ===
              complex(Inf)
    end

    @testset "the sentinel" begin
        LFTs = LinearFractionalTransformations

        @test ∞ isa RiemannInfinity
        @test LFTs.∞ === RiemannInfinity()
        @test sprint(show, ∞) == "∞"
        @test isinf(∞)
        @test LFTs.is_infinite_point(∞)
        @test LFTs.is_infinite_point(-Inf)
        @test LFTs.is_infinite_point(complex(0, Inf))
        @test !LFTs.is_infinite_point(1 // 2)

        # ∞ is accepted wherever Inf is, but never as a coefficient
        h = LFT(∞, 1.0 + 0im, 2.0 + 0im)
        @test h(∞) == 0
        @test h(1.0 + 0im) == 1
        @test isinf(h(2.0 + 0im))
        @test_throws ErrorException LFT(∞, 0, 0, 1)
    end

    @testset "a field with no isinf" begin
        at(x) = ToyExact(x)
        m = LFT(at(1), at(2), at(3), at(4))
        @test m(at(0)) == at(1 // 2)
        @test m(∞) == at(1 // 3)

        poled = LFT(at(0), at(1), at(1), at(0))          # z -> 1/z
        @test poled(at(2)) == at(1 // 2)
        @test poled(at(0)) === ∞                         # the fallback
        @test LinearFractionalTransformations.infinity(at(0)) === ∞
        @test map(poled, [at(0), at(2)]) isa
              Vector{Union{ToyExact, RiemannInfinity}}
    end

    @testset "hash" begin
        @test hash(LFT(1.0, 0.0, 0.0, 1.0)) isa UInt
        @test hash(LFT(0.0, 1.0, 1.0, 0.0)) ==
              hash(LFT(0 // 1, 1 // 1, 1 // 1, 0 // 1))

        at(x) = ToyExact(x)
        poled = LFT(at(0), at(1), at(1), at(0))
        scaled = LFT(at(0), at(2), at(2), at(0))
        @test hash(poled) isa UInt
        @test hash(poled) == hash(scaled)
    end
end
