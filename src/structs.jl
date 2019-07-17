mutable struct SourceTerm
    λ::Float64
    μ::Float64
end

mutable struct DTRWSolver
    a::Float64
    b::Float64
    x::Array{Float64,1}
    mid_price::Float64
end

mutable struct ReactionDiffusionPricePath
    T::Int64
    τ::Int64
    initial_mid_price::Float64
    n_spatial_points::Int64
    β::Float64
    sample_std::Float64
    σ::Float64
    D::Float64
    ν::Float64
    α::Float64
    source_term::SourceTerm
end
ReactionDiffusionPricePath(dict)=ReactionDiffusionPricePath(
    dict["T"], dict["τ"], dict["initial_mid_price"], dict["n_spatial_points"],
    dict["boltz_const"], dict["sample_std"], dict["σ"], dict["D"], dict["ν"],
    dict["α"], SourceTerm(dict["λ"], dict["μ"]))

ReactionDiffusionPricePath(T::Int64=100, τ::Int64=10,
    initial_mid_price::Real=100.0, n_spatial_points::Int64=101,
    boltz_const::Real=1.0, sample_std::Real=4.0, σ::Real=0.001, D::Real=5.0,
    ν::Real=0.001, α::Real=1.0, λ::Real=1.0, μ::Real=0.5) =
    ReactionDiffusionPricePath(T, τ, initial_mid_price,
    n_spatial_points, boltz_const, sample_std, σ, D, ν, α,
    SourceTerm(λ, μ))
