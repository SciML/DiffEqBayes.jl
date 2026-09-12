using DiffEqBayes, BenchmarkTools
using OrdinaryDiffEqTsit5, Distributions, Turing
using StableRNGs

const SUITE = BenchmarkGroup()

rng = StableRNG(1234)

function lotka_volterra!(du, u, p, t)
    x, y = u
    a = p[1]
    du[1] = a * x - x * y
    du[2] = -3y + x * y
    return nothing
end

u0 = [1.0, 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lotka_volterra!, u0, tspan, [1.5])
sol = solve(prob, Tsit5())
t = collect(range(1; stop = 10, length = 10))
data = reduce(hcat, [sol(ti) .+ 0.01 .* randn(rng, 2) for ti in t])
priors = [Normal(1.5, 0.5)]

# =============================================================================
# Stan code generation + data plumbing (no cmdstan needed)
# =============================================================================

SUITE["stan"] = BenchmarkGroup()

SUITE["stan"]["stan_string"] = @benchmarkable DiffEqBayes.stan_string($(Normal(0.0, 1.0)))
SUITE["stan"]["stan_ode_data"] = @benchmarkable StanODEData()

# =============================================================================
# Turing inference — the representative workload: NUTS sampling over ODE
# solutions. Bounded by num_samples.
# =============================================================================

SUITE["inference"] = BenchmarkGroup()

SUITE["inference"]["turing_nuts"] = @benchmarkable turing_inference(
    $prob, Tsit5(), $t, $data, $priors;
    sample_args = (num_samples = 60,), progress = false
) seconds = 900
