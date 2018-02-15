get_plot_object(bayesian_result::StanModel) = Mamba.plot(bayesian_result.chain_results)
get_plot_object(bayesian_result::Turing.Chain) = Mamba.plot(bayesian_result)

function plot_chain(bayesian_result::StanModel,filename=nothing)
	p = get_plot_object(bayesian_result)
	if filename == nothing
		return Mamba.draw(p)
	else
		return Mamba.draw(p, filename=filename)
	end
end

function plot_chain(bayesian_result::Turing.Chain,filename=nothing)
	p = Mamba.plot(bayesian_result)
	if filename == nothing
		return Mamba.draw(p)
	else
		return Mamba.draw(p, filename=filename)
	end
end