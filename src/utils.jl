get_plot_object(bayesian_result::StanModel) = Mamba.plot(bayesian_result.chain_results)

function plot_chain(bayesian_result::StanModel,filename=nothing)
	p = get_plot_object(bayesian_result)
	if filename == nothing
		return Mamba.draw(p)
	else
		return Mamba.draw(p, filename=filename)
	end
end

function plot_chain(bayesian_result::T,filename=nothing) where {T<:MCMCChains.Chains}
	StatsPlots.plot(bayesian_result)
end