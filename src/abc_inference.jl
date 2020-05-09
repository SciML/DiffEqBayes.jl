function createabcfunction(prob, t, distancefunction, alg; save_idxs = nothing, sample_u0 = false, kwargs...)
    function simfunc(params, constants, data)
        local u0
        if sample_u0
            u0 = save_idxs === nothing ? params[1:length(prob.u0)] : params[1:length(save_idxs)]
            if length(u0) < length(prob.u0)
                for i in length(u0):length(prob.u0)
                    push!(u0,prob.u0[i])
                end
            end
        else
            u0 = prob.u0
        end
        sol = concrete_solve(STANDARD_PROB_GENERATOR(prob, params), alg, u0; saveat = t, save_idxs = save_idxs, kwargs...)
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
                       num_samples = 500, maxiterations = 10^5, save_idxs = nothing, sample_u0 = false, kwargs...)

    abcsetup = ABCalgorithm(createabcfunction(prob, t, distancefunction, alg; save_idxs = save_idxs, sample_u0 = sample_u0, kwargs...),
                            length(priors),
                            ϵ,
                            ApproxBayes.Prior(priors);
                            nparticles = num_samples,
                            maxiterations = maxiterations
                            )

    abcresult = runabc(abcsetup, data, progress = progress)
    return abcresult
end
