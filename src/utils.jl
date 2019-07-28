function plot_chain(bayesian_result::T,filename=nothing) where {T<:MCMCChains.Chains}
	StatsPlots.plot(bayesian_result)
end
