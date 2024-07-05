#--- general ---#
x0 = 0
x1 = 1
nx = 5

KC.linspace(x0, x1, nx)
KC.heaviside(x0)
KC.fortsign(x0, x1)

KC.convergence_order(1e-2, 1e-3)
KC.L1_error(rand(3), rand(3), 1e-2)
KC.L2_error(rand(3), rand(3), 1e-2)
KC.L∞_error(rand(3), rand(3), 1e-2)

KC.lgwt(12, -1, 1)

#--- polylog ---#
KC.polylog_duplication(1, 1)
KC.polylog_duplication(1 + 1im, 1)

# check throws errors
@testset "throws errors" begin
    @test_throws MethodError KC.polylog(1, "1.0")  # should work for all numbers
    @test_throws MethodError KC.polylog("1.0", 1)  # should work for all numbers
    @test_throws MethodError KC.polylog(1)
end

# check output types
@testset "output types" begin
    @test typeof(KC.polylog(1, 0)) == Float64
    @test typeof(KC.polylog(complex(1), 0)) == Complex{Float64}
end

# check keywords work
@testset "output types" begin
    @test KC.polylog(
        -1,
        0;
        level = 1,
        accuracy = 1.0e-5,
        min_iterations = 1,
        max_iterations = 100,
    ) ≈ 0.0
end

# check subsidary routine errors, just to make coverage cleaner
@testset "unexported function errors" begin
    @test_throws DomainError KC.polylog_series_1(1.0, 2.0)
    @test_throws DomainError KC.polylog_series_1(1.0, 0.75)
    @test_throws DomainError KC.polylog_series_2(1.0, 0.0)
    @test_throws DomainError KC.polylog_series_3(1.0, 0.0)
    @test_throws DomainError KC.polylog_series_3(-1.0, 0.6)
    @test_throws DomainError KC.g_crandall(-1)
    @test_throws DomainError KC.g_crandall(15)
    τ = 0.00001
    @test_throws DomainError KC.Q_closed(0, τ, 0; n_terms = 0)
    @test_throws DomainError KC.Q_closed(0, τ, 0; n_terms = 4)
    @test_throws DomainError KC.Q(0, τ, 0; n_terms = 0)
    @test_throws DomainError KC.Q(0, τ, 0; n_terms = 8)
end

KC.c_closed(0, 0, 0.1)
KC.c_closed(0, 1, 0.1)
KC.c_closed(0, 2, 0.1)
KC.Q(0, 0.1, 0; n_terms = 4)
KC.c_crandall(0, 0, 0)
KC.b_crandall(0, 0, 0)
KC.f_crandall(0, 0)
KC.harmonic(1)
KC.harmonic(0, 1.0)
KC.harmonic(1, 1.0)
KC.harmonic(2, 3.0)
KC.harmonic(2, 3)
KC.stieltjes(2)
