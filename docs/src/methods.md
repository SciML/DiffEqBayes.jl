# Bayesian Methods

The following methods require DiffEqBayes.jl. The backend packages are loaded as
dependencies of DiffEqBayes, but Stan also requires a local CmdStan installation.

```julia
using Pkg
Pkg.add("DiffEqBayes")
using DiffEqBayes
```

The API pages below describe the current signatures, keyword defaults, return values,
and backend-specific constraints. The older ApproxBayes interface is not part of the
current package: its implementation remains disabled and is intentionally omitted here.

## Stan

[`stan_inference`](@ref) uses
[StanSample.jl](https://stanjulia.github.io/StanSample.jl/stable/) for Bayesian
inference. See the [Stan installation guide](https://stanjulia.github.io/StanSample.jl/stable/INSTALLATION/)
before using it.

## Turing

[`turing_inference`](@ref) uses
[Turing.jl](https://turinglang.org/stable/) and returns the chain produced by
`Turing.sample`. Parameter draws are indexed by the `VarName` values supplied through
`syms`.

## DynamicHMC

[`dynamichmc_inference`](@ref) uses
[DynamicHMC.jl](https://www.tamaspapp.eu/DynamicHMC.jl/stable/) and returns the
backend result with a `posterior` field containing transformed parameter values.

## Docstrings

```@docs
DiffEqBayes
StanODEData
StanResult
stan_inference
turing_inference
dynamichmc_inference
```
