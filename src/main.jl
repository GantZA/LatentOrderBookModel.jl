using Plots

include("reaction_diffusion_path.jl")
include("preprocess_bar_data.jl")


bar_data = process_bar_data("Original_Price_Bars_2300.csv")

sim_price_path = reaction_diffusion_path(bar_data, 0.0003, 5, 0.001,
    source_term, [1.0, 0.01])


plot(1:2300, sim_price_path)
plot!(1:2300, bar_data.price)
