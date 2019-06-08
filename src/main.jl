using Plots

include("reaction_diffusion_path.jl")
include("preprocess_bar_data.jl")


bar_data = process_bar_data("Original_Price_Bars_2300.csv")
seed = 40
Random.seed!(seed)

mutable struct RDP_Params
    T::Int64
    τ::Int64
    initial_mid_price::Float64
    n_spatial_points::Int64
    boltz_const::Float64
    sample_std::Float64
    σ::Float64
    D::Float64
    η::Float64
end

mutable struct ST_Params
    λ::Float64
    μ::Float64
end

rdp_params = RDP_Params()
st_params = ST_Params()

sim_price_path = reaction_diffusion_path(
    rd_params,
    st_params)


plot(1:2300, sim_price_path)
plot!(1:2300, bar_data.price)


mutable struct test111
    a::Float64
    b::Float64
end

a = test111(4,4212)
typeof(a)
