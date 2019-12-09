# LOB
mutable struct ObjectiveSurface
    param1_name::String
    param1_values::Array{Any, 1}
    param2_name::String
    param2_values::Array{Any, 1}
    surface_points::Int64
    replications::Int64
    parameters::Dict{String, Any}
end
function ObjectiveSurface(param1_name::String, param1_values::Array{Float64, 1},
    param2_name::String, param2_values::Array{Float64, 1},
    surface_points::Int64, parameters::Dict)
    return ObjectiveSurface(param1_name, param1_values, param2_name,
    param2_values, surface_points, parameters["num_paths"], parameters)
end


function (os::ObjectiveSurface)(seed, para=false)
    iterations = os.surface_points^2

    seeds = Int.(rand(MersenneTwister(seed), UInt32, iterations))
    if para==true
        sample_price_paths = SharedArray{Float64,2}((os.parameters["T"]+1, iterations*os.replications))
        @sync @distributed for i in 1:iterations
            os.parameters[os.param1_name] = os.param1_values[i]
            os.parameters[os.param2_name] = os.param2_values[i]
            rdpp = ReactionDiffusionPricePaths(os.parameters)
            index_start = (i-1) * os.replications + 1
            index_end = i * os.replications
            _, _, sample_price_paths[:, index_start:index_end], _, _ = rdpp(seeds[i])
        end
        return sample_price_paths
    else
        sample_price_paths = Array{Float64,2}(undef, os.parameters["T"]+1, iterations*os.replications)
        for i in 1:iterations
            os.parameters[os.param1_name] = os.param1_values[i]
            os.parameters[os.param2_name] = os.param2_values[i]
            rdpp = ReactionDiffusionPricePaths(os.parameters)
            index_start = (i-1) * os.replications + 1
            index_end = i * os.replications
            _, _, sample_price_paths[:, index_start:index_end], _, _ = rdpp(seeds[i])
        end
        return sample_price_paths
    end
end
