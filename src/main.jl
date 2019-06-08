using Plots
using Random

include("reaction_diffusion_path.jl")
include("preprocess_bar_data.jl")


bar_data = process_bar_data("data/Original_Price_Bars_2300.csv")


# mutable struct RDP_Params
#     T::Int64
#     τ::Int64
#     initial_mid_price::Float64
#     n_spatial_points::Int64
#     boltz_const::Float64
#     sample_std::Float64
#     σ::Float64
#     D::Float64
#     η::Float64
# end

# mutable struct ST_Params
#     λ::Float64
#     μ::Float64
# end
bar_data.price
using Statistics
std(bar_data.price)
rdp_params = RDP_Params(
    2300, 10, 238.745,
    501, 2, 7.415, 0.001,
    5.0, 0.001)
st_params = ST_Params(1.0, 0.5)

seed = 445784574
Random.seed!(seed)

sim_price_path = reaction_diffusion_path(
    rdp_params,
    st_params)


plot(1:2300, sim_price_path)
plot!(1:2300, bar_data.price)
