function (os::ObjectiveSurface)(seed)
    iterations = os.surface_points^2 * os.replications

    seeds = Int.(rand(MersenneTwister(seed), UInt32, os.replications, iterations))
    price_paths = Array{Float64,2}(undef, os.params["T"], iterations)
    for i in 1:iterations
        os.params[os.param1_name] = os.param1_values[i]
        os.params[os.param2_name] = os.param2_values[i]
        rdpp = ReactionDiffusionPricePath(os.params)
        price_paths[:,i] = rdpp(seeds[i])
    end
    return price_paths
end
