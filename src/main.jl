using ArgParse

include("reaction_diffusion_path.jl")
# include("preprocess_bar_data.jl")


# bar_data = process_bar_data("data/Original_Price_Bars_2300.csv")


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

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "seed"
            help = "Seed for randomness"
            required = true
            arg_type = Int
        "T"
            help = "Number of time periods"
            required = true
            arg_type = Int
        "τ"
            required = true
            arg_type = Int
        "initial_mid_price"
            required = true
            arg_type = Float64
        "n_spatial_points"
            required = true
            arg_type = Int
        "boltz_const"
            required = true
            arg_type = Float64
        "sample_std"
            required = true
            arg_type = Float64
        "σ"
            required = true
            arg_type = Float64
        "D"
            required = true
            arg_type = Float64
        "η"
            required = true
            arg_type = Float64
        "λ"
            required = true
            arg_type = Float64
        "μ"
            required = true
            arg_type = Float64
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    rdp_params = RDP_Params(
        parse_args["seed"],
        parsed_args["T"], parsed_args["τ"], parsed_args["initial_mid_price"],
        parsed_args["n_spatial_points"], parsed_args["boltz_const"],
        parsed_args["sample_std"], parsed_args["σ"], parsed_args["D"],
        parsed_args["η"])
    st_params = ST_Params(parsed_args["λ"], parsed_args["μ"])
    print(reaction_diffusion_path(
        rdp_params,
        st_params))
end


# rdp_params = RDP_Params(
#     2300, 10, 238.745, 501, 2, 7.415, 0.001, 5.0, 0.001, 1.0, 0.5)
# st_params = ST_Params(1.0, 0.5)
#
# seed = 445784574
# Random.seed!(seed)
#
# sim_price_path = reaction_diffusion_path(
#     rdp_params,
#     st_params)
#
#
# plot(1:2300, sim_price_path)
# plot!(1:2300, bar_data.price)
