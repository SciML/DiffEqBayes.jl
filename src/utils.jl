using Mamba
using DiffEqBayes

function plot_chain(bayesian_result::StanModel,filename=nothing)
	p = Mamba.plot(bayesian_result.chain_results)
	if filename == nothing
		return Mamba.draw(p)
	else
		return Mamba.draw(p, filename=filename)
	end
end