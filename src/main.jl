using LatentOrderBookModel

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

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    main()
end
