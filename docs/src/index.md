# DiffEqBayes.jl: Bayesian Parameter Estimation for Differential Equations

This repository is a set of extension functionality for estimating the parameters
of differential equations using Bayesian methods. It allows the choice of using
[CmdStan.jl](https://stanjulia.github.io/CmdStan.jl/stable/), [Turing.jl](https://turing.ml/stable/docs/using-turing/), [DynamicHMC.jl](https://www.tamaspapp.eu/DynamicHMC.jl/stable/) and
[ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl) to perform a
Bayesian estimation of a differential equation problem specified via the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) interface.

## Installation

To install DiffEqBayes.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("DiffEqBayes")
```

## Contributing

- Please refer to the
  [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to SciML.
- See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
- There are a few community forums:
    - The #diffeq-bridged and #sciml-bridged channels in the
      [Julia Slack](https://julialang.org/slack/)
    - The #diffeq-bridged and #sciml-bridged channels in the
      [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
    - On the [Julia Discourse forums](https://discourse.julialang.org)
    - See also [SciML Community page](https://sciml.ai/community/)
