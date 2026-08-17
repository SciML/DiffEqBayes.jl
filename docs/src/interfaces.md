# Public Interfaces

The inference entry points accept SciML problems through the common problem interface.
The rules below are the contract used by the implementations and by the generic tests
in `test/public_interface.jl`.

## Problem interface

`stan_inference`, `turing_inference`, and `dynamichmc_inference` accept an
`SciMLBase.AbstractSciMLProblem`. For Turing and DynamicHMC, the problem must support
`solve(problem, algorithm; ...)`. Stan instead translates the problem using its standard
SciML fields `u0`, `p`, and, when required by the algorithm, `tspan`.

The inference code may replace `u0` and `p` for every sampled parameter vector. The
problem therefore must support the normal SciML parameter contract: its parameter
container must accept the sampled values, and `SciMLStructures.Tunable` must be
available when a structured parameter container is used. Plain arrays and ordinary
SciML problem types are the simplest supported form.

## Observations

For time-dependent problems, `t` contains the observation times and `data[:, i]`
contains the observation at `t[i]`. The number of columns in `data` must match the
number of entries in `t`, and the number of rows must match the observed state
components. The same state ordering must be used by `save_idxs` when only part of the
state is observed.

For `turing_inference`, `t = nothing` is supported for a one-dimensional result such as
a steady-state problem. In that case `data` is passed as one observation vector. The
DynamicHMC and Stan interfaces require an explicit time grid.

## Priors and parameters

`priors` has one entry for every sampled parameter. With `sample_u0 = true`, the
selected initial-condition entries precede the model parameters in the sampled vector.
The names in `syms` must have the same length and identify those sampled entries.

For DynamicHMC, `parameter_transformations` must map an unconstrained real vector to
the valid parameter container, and `mcmc_kwargs.initialization.q` must contain one
entry for every sampled parameter and noise scale. `σ_priors` has one entry for each
observed component.

## Backend rules

- `turing_inference` forwards `solve_kwargs` to SciML `solve`, and forwards sampling
  options through `sample_args` and `sample_kwargs`. The `likelihood` callable receives
  `(u, p, t, σ)` and must return the distribution used for the observation.
- `dynamichmc_inference` uses `solve_kwargs` for SciML `solve` and
  `mcmc_kwargs` for DynamicHMC. Its result contains the backend result plus a
  `posterior` collection of transformed parameter values.
- `stan_inference` supports Stan's `:adams`, `:rk45`, and `:bdf` algorithm names when
  generating a model. `solve_kwargs`, `sample_kwargs`, and `output_format` are passed
  to the corresponding Stan operations. It requires a CmdStan installation.

## Extension rules

The supported extension point is a new `SciMLBase.AbstractSciMLProblem` together with
the ordinary SciML `solve` interface. Users should call the exported inference
functions and should not construct or extend `DynamicHMCPosterior`; that type is an
implementation detail used to adapt the DynamicHMC backend. Backend-specific solver
and sampler options should be passed through the documented keyword containers rather
than by reaching into DiffEqBayes internals. A new problem type must therefore be
usable through the SciMLBase problem and `solve` interfaces before it can be used with
these entry points; adding a DiffEqBayes-specific method for an internal helper is not
a supported extension.

The generic contract test in `test/public_interface.jl` uses only the exported
inference function, a public `SciMLBase.ODEProblem`, and the public `solve`-compatible
solver interface. It deliberately does not call `DynamicHMCPosterior`, internal
likelihood closures, or generated backend models.

## Developer Interface

`DynamicHMCPosterior` is documented for maintainers and backend developers who need to
understand the adapter used by `dynamichmc_inference`. It is not user-facing API:
applications should call [`dynamichmc_inference`](@ref) and should not construct,
subtype, or extend this type.

```@docs
DiffEqBayes.DynamicHMCPosterior
```
