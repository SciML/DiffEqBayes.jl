# DiffEqBayes.jl

This repository is a set of extension functionality for estimating the parameters
of differential equations using Bayesian methods. It allows the choice of using
[CmdStan.jl](https://github.com/StanJulia/CmdStan.jl), [Turing.jl](https://github.com/TuringLang/Turing.jl), [DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl) and
[ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl) to perform a
Bayesian estimation of a differential equation problem specified via the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) interface.

## Installation

For the Bayesian methods, you must install DiffEqBayes.jl:

```julia
]add DiffEqBayes
using DiffEqBayes
```
