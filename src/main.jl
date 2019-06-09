using ArgParse
using LatentOrderBookModel

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "SEED"
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
        parsed_args["SEED"],
        parsed_args["T"], parsed_args["τ"], parsed_args["initial_mid_price"],
        parsed_args["n_spatial_points"], parsed_args["boltz_const"],
        parsed_args["sample_std"], parsed_args["σ"], parsed_args["D"],
        parsed_args["η"])
    st_params = ST_Params(parsed_args["λ"], parsed_args["μ"])
    print(reaction_diffusion_path(
        rdp_params,
        st_params))
end

main()
