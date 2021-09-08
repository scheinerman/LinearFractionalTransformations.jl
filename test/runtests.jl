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

end


# check out LFTQ 
@testset "LFTQ" begin
    for k = 1:5
        Q = [1 0; 0 -1]
        while det(Q) < 0.5
            Q, R = qr(randn(3, 3))
        end
        F = LFTQ(Q)
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
