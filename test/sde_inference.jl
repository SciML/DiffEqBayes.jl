using DiffEqBayes
using DiffEqBayes: _em_loglikelihood
using Distributions, ForwardDiff, LinearAlgebra, Random
using SciMLBase: SDEProblem
using Test
using Turing

# Central finite differences, used as an AD-independent gradient reference.
function central_gradient(f, x; h = 1.0e-6)
    return map(eachindex(x)) do j
        step = zeros(length(x))
        step[j] = h
        return (f(x .+ step) - f(x .- step)) / (2h)
    end
end

@testset "Scalar SDE with state-dependent diagonal noise" begin
    drift!(du, u, p, t) = (du[1] = p[1] * (p[2] - u[1]))
    diffusion!(du, u, p, t) = (du[1] = p[3] * sqrt(abs(u[1]) + 1))
    drift(u, p, t) = [p[1] * (p[2] - u[1])]
    diffusion(u, p, t) = [p[3] * sqrt(abs(u[1]) + 1)]
    scalar_drift(u, p, t) = p[1] * (p[2] - u)
    scalar_diffusion(u, p, t) = p[3] * sqrt(abs(u) + 1)

    t = [0.0, 0.1, 0.3, 0.6]
    data = [0.5 0.7 0.4 0.9]
    p = [1.2, 0.8, 0.6]
    inplace_prob = SDEProblem(drift!, diffusion!, [0.5], (0.0, 0.6), p)
    outofplace_prob = SDEProblem(drift, diffusion, [0.5], (0.0, 0.6), p)
    scalar_prob = SDEProblem(scalar_drift, scalar_diffusion, 0.5, (0.0, 0.6), p)

    # Independent evaluation of the Euler-Maruyama transition densities.
    function reference(p, sigma)
        return sum(1:(length(t) - 1)) do i
            dt = t[i + 1] - t[i]
            x = data[1, i]
            mean_next = x + p[1] * (p[2] - x) * dt
            variance = abs2(p[3]) * (abs(x) + 1) * dt + abs2(sigma)
            return logpdf(Normal(mean_next, sqrt(variance)), data[1, i + 1])
        end
    end

    # Reference values computed separately with SciPy.
    @test reference(p, 0.0) ≈ -0.399584716884109
    @test reference(p, 0.2) ≈ -0.6362096808210967

    for prob in (inplace_prob, outofplace_prob, scalar_prob), sigma in (0.0, 0.2)
        @test _em_loglikelihood(prob, p, t, data, sigma) ≈ reference(p, sigma)
        # Repeated evaluations target one fixed density: nothing is simulated.
        @test _em_loglikelihood(prob, p, t, data, sigma) ==
            _em_loglikelihood(prob, p, t, data, sigma)
        loglikelihood = q -> _em_loglikelihood(prob, q, t, data, sigma)
        @test ForwardDiff.gradient(loglikelihood, p) ≈
            central_gradient(loglikelihood, p) rtol = 1.0e-6
    end

    # The diffusion parameter enters the likelihood.
    @test _em_loglikelihood(inplace_prob, [1.2, 0.8, 0.3], t, data, 0.0) ≈
        -1.611026461852265

    # A noise-free component has no transition density unless `sigma` regularizes it.
    @test _em_loglikelihood(inplace_prob, [1.2, 0.8, 0.0], t, data, 0.0) == -Inf
    @test isfinite(_em_loglikelihood(inplace_prob, [1.2, 0.8, 0.0], t, data, 0.2))
end

@testset "Non-diagonal noise uses the full transition covariance" begin
    # Two states driven by one shared Wiener process: perfectly correlated increments.
    no_drift!(du, u, p, t) = fill!(du, 0)
    shared_noise!(du, u, p, t) = fill!(du, p[1])
    shared_prob = SDEProblem(
        no_drift!, shared_noise!, zeros(2), (0.0, 1.0), [1.0];
        noise_rate_prototype = zeros(2, 1)
    )
    shared_data = [0.0 1.0; 0.0 -1.0]
    shared = _em_loglikelihood(shared_prob, [1.0], [0.0, 1.0], shared_data, 0.1)
    @test shared ≈ logpdf(MvNormal(zeros(2), [1.01 1.0; 1.0 1.01]), [1.0, -1.0])
    @test shared ≈ -99.88435933445069
    # Treating the two components as independently driven would give about -2.84.
    @test shared < -50

    # Two states and three Wiener processes; `noise!` writes only the nonzero entries.
    function drift!(du, u, p, t)
        du[1] = -p[1] * u[1]
        du[2] = u[1] - u[2]
        return nothing
    end
    function noise!(du, u, p, t)
        du[1, 1] = p[2]
        du[1, 3] = p[3] * u[2]
        du[2, 2] = p[4]
        du[2, 3] = p[3]
        return nothing
    end
    drift(u, p, t) = [-p[1] * u[1], u[1] - u[2]]
    noise(u, p, t) = [p[2] 0 p[3] * u[2]; 0 p[4] p[3]]

    t = [0.0, 0.2, 0.5, 0.6]
    data = [1.0 1.3 0.8 1.1; 0.5 0.2 0.6 0.4]
    p = [0.7, 0.4, 0.3, 0.5]
    inplace_prob = SDEProblem(
        drift!, noise!, [1.0, 0.5], (0.0, 0.6), p; noise_rate_prototype = zeros(2, 3)
    )
    outofplace_prob = SDEProblem(
        drift, noise, [1.0, 0.5], (0.0, 0.6), p; noise_rate_prototype = zeros(2, 3)
    )

    function reference(p, sigma)
        return sum(1:(length(t) - 1)) do i
            dt = t[i + 1] - t[i]
            u = data[:, i]
            G = noise(u, p, t[i])
            covariance = Matrix(Symmetric(G * G' * dt + Diagonal(abs2.(sigma))))
            mean_next = u + drift(u, p, t[i]) * dt
            return logpdf(MvNormal(mean_next, covariance), data[:, i + 1])
        end
    end

    # Reference values computed separately with SciPy.
    @test reference(p, [0.0, 0.0]) ≈ -6.2130645449992175
    @test reference(p, [0.1, 0.2]) ≈ -3.7725678551841892

    for prob in (inplace_prob, outofplace_prob)
        @test _em_loglikelihood(prob, p, t, data, 0.0) ≈ reference(p, [0.0, 0.0])
        @test _em_loglikelihood(prob, p, t, data, [0.1, 0.2]) ≈ reference(p, [0.1, 0.2])
        loglikelihood = q -> _em_loglikelihood(prob, q, t, data, [0.1, 0.2])
        @test ForwardDiff.gradient(loglikelihood, p) ≈
            central_gradient(loglikelihood, p) rtol = 1.0e-6
    end
end

@testset "Invalid inputs are rejected before sampling" begin
    drift!(du, u, p, t) = (du .= -p[1] .* u)
    diffusion!(du, u, p, t) = fill!(du, p[2])
    prob = SDEProblem(drift!, diffusion!, zeros(2), (0.0, 1.0), [1.0, 0.5])
    priors = [Uniform(0, 5), Uniform(0, 2)]
    t = [0.0, 0.5, 1.0]
    data = [0.0 0.1 0.2; 0.0 -0.1 0.3]

    @test_throws ArgumentError sde_inference(prob, [0.0, 0.5, 0.5], data, priors)
    @test_throws ArgumentError sde_inference(prob, [0.0], data[:, 1:1], priors)
    @test_throws DimensionMismatch sde_inference(prob, t, data[:, 1:2], priors)
    @test_throws DimensionMismatch sde_inference(prob, t, data[1:1, :], priors)
    @test_throws ArgumentError sde_inference(prob, t, data, priors; sigma = -0.1)
    @test_throws DimensionMismatch sde_inference(prob, t, data, priors; sigma = [0.1])

    matrix_prob = SDEProblem(drift!, diffusion!, zeros(2, 1), (0.0, 1.0), [1.0, 0.5])
    @test_throws ArgumentError sde_inference(matrix_prob, t, data, priors)

    # Fewer Wiener processes than states makes the covariance singular without `sigma`.
    shared_prob = SDEProblem(
        drift!, diffusion!, zeros(2), (0.0, 1.0), [1.0, 0.5];
        noise_rate_prototype = zeros(2, 1)
    )
    @test_throws ArgumentError sde_inference(shared_prob, t, data, priors)
end

@testset "Sampling through the exported interface" begin
    # Out-of-place scalar problem with observations passed as a vector.
    scalar_prob = SDEProblem(
        (u, p, t) -> p[1] * (p[2] - u), (u, p, t) -> p[3], 0.0, (0.0, 1.0),
        [1.0, 0.0, 0.5]
    )
    chain = sde_inference(
        scalar_prob, [0.0, 0.5, 1.0], [0.0, 0.2, 0.1],
        [Uniform(0, 5), Normal(0, 1), Uniform(0.05, 2)];
        syms = [:rate, :level, :noise],
        sample_args = (sampler = Turing.MH(), num_samples = 2)
    )
    @test size(chain, 1) > 0

    # Ornstein-Uhlenbeck process simulated from its exact transition density, so the
    # data do not depend on the discretization being tested.
    rate, level, noise = 1.0, 2.0, 0.5
    dt = 0.1
    t = collect(0.0:dt:99.9)
    rng = Random.Xoshiro(2026)
    data = Matrix{Float64}(undef, 1, length(t))
    data[1, 1] = level
    decay = exp(-rate * dt)
    stationary_step = noise * sqrt((1 - abs2(decay)) / (2 * rate))
    for i in 1:(length(t) - 1)
        data[1, i + 1] = level + (data[1, i] - level) * decay +
            stationary_step * randn(rng)
    end

    drift!(du, u, p, t) = (du[1] = p[1] * (p[2] - u[1]))
    diffusion!(du, u, p, t) = (du[1] = p[3])
    prob = SDEProblem(drift!, diffusion!, [level], (0.0, 99.9), [rate, level, noise])
    priors = [Uniform(0, 5), Normal(0, 5), Uniform(0.05, 2)]

    Random.seed!(1234)
    chain = sde_inference(
        prob, t, data, priors;
        syms = [:rate, :level, :noise], sample_args = (num_samples = 500,)
    )
    @test mean(chain[@varname(rate)]) ≈ rate atol = 0.6
    @test mean(chain[@varname(level)]) ≈ level atol = 0.3
    # The diffusion parameter is identified by the data, not left at its prior mean.
    @test mean(chain[@varname(noise)]) ≈ noise atol = 0.1
end
