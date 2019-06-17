module LatentOrderBookModel

include("reaction_diffusion_path.jl")
include("parse_params.jl")


reaction_diffusion_path(SEED::Int64=1, T::Int64=100, τ::Int64=10,
    initial_mid_price::Real=100.0, n_spatial_points::Int64=101, boltz_const::Real=2.0,
    sample_std::Real=4.0, σ::Real=0.001, D::Real=5.0, η::Real=0.001, λ::Real=1.0,
    μ::Real=0.5) =
    reaction_diffusion_path(RDP_Params(SEED, T, τ, initial_mid_price,
    n_spatial_points, boltz_const, sample_std, σ, D, η), ST_Params(λ, μ))

export reaction_diffusion_path, RDP_Params, ST_Params, parse_commandline

end # module
