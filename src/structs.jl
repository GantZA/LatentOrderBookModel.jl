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
    D::Float64
    ν::Float64
    α::Float64
    source_term::SourceTerm
end
ReactionDiffusionPricePath(dict)=ReactionDiffusionPricePath(
    dict["T"], dict["τ"], dict["initial_mid_price"], dict["n_spatial_points"],
    dict["boltz_const"], dict["sample_std"], dict["D"], dict["ν"],
    dict["α"], SourceTerm(dict["λ"], dict["μ"]))

ReactionDiffusionPricePath(;T::Int64=100, τ::Int64=10,
    initial_mid_price::Real=100.0, n_spatial_points::Int64=101,
    boltz_const::Real=1.0, sample_std::Real=4.0, D::Real=5.0,
    ν::Real=0.001, α::Real=1.0, λ::Real=1.0, μ::Real=0.5) =
    ReactionDiffusionPricePath(T, τ, initial_mid_price,
    n_spatial_points, boltz_const, sample_std, D, ν, α,
    SourceTerm(λ, μ))


mutable struct ObjectiveSurface
    param1_name::String
    param1_values::Array{Any,1}
    param2_name::String
    param2_values::Array{Any,1}
    surface_points::Int64
    replications::Int64
    params::Dict{String, Any}
end
