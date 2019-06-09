module LatentOrderBookModel

include("reaction_diffusion_path.jl")

reaction_diffusion_path(SEED::Int64, T::Int64, τ::Int64,
    initial_mid_price::Float64, n_spatial_points::Int64, boltz_const::Int64,
    sample_std::Float64, σ::Float64, D::Float64, η::Float64, λ::Float64,
    μ::Float64) =
    reaction_diffusion_path(RDP_Params(SEED, T, τ, initial_mid_price,
    n_spatial_points, boltz_const, sample_std, σ, D, η), ST_Params(λ, μ))

export reaction_diffusion_path

end # module
