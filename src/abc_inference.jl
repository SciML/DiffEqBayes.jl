function createabcfunction(prob, t, distancefunction, alg; save_idxs = nothing, kwargs...)
    function simfunc(params, constants, data)
        sol = concrete_solve(STANDARD_PROB_GENERATOR(prob, params), alg; saveat = t, save_idxs = save_idxs, kwargs...)
        if size(sol, 2) < length(t)
            return Inf,nothing
        else
            simdata = convert(Array, sol)
            return distancefunction(data, simdata), nothing
        end
    end
end

function abc_inference(prob::DiffEqBase.DEProblem, alg, t, data, priors; ϵ=0.001,
                       distancefunction = euclidean, ABCalgorithm = ABCSMC, progress = false,
                       num_samples = 500, maxiterations = 10^5, save_idxs = nothing, kwargs...)

    abcsetup = ABCalgorithm(createabcfunction(prob, t, distancefunction, alg; save_idxs = save_idxs, kwargs...),
                            length(priors),
                            ϵ,
                            ApproxBayes.Prior(priors);
                            nparticles = num_samples,
                            maxiterations = maxiterations
                            )

    abcresult = runabc(abcsetup, data, progress = progress)
    return abcresult
end
