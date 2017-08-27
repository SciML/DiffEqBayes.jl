# DiffEqBayes.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/DiffEqBayes.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DiffEqBayes.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/DiffEqBayes.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/DiffEqBayes.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/DiffEqBayes.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/DiffEqBayes.jl?branch=master)

This repository is a set of extension functionality for estimating the parameters of differential equations using Bayesian methods. It allows the choice of using [Stan.jl](https://github.com/goedman/Stan.jl) and [Turing.jl](https://github.com/yebai/Turing.jl) to perform a Bayesian estimation of a differential equation problem specified via the [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) interface. 
