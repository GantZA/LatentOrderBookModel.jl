using LatentOrderBookModel

function main(output="stdout")
    parsed_args = parse_commandline()
    rdpp = ReactionDiffusionPricePaths(parsed_args["T"], parsed_args["τ"],
        parsed_args["initial_mid_price"], parsed_args["n_spatial_points"],
        parsed_args["boltz_const"], parsed_args["sample_std"],
        parsed_args["D"], parsed_args["nu"], parsed_args["σ"],
        SourceTerm(parsed_args["λ"], parsed_args["μ"]))
    if output=="stdout"
        print(rdpp(parsed_args["SEED"]))
    else
        return rdpp(parsed_args["SEED"])
    end


end

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    main()
end
