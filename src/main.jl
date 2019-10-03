using LatentOrderBookModel

function main(output="stdout")
    parsed_args = parse_commandline()
    rdpp = ReactionDiffusionPricePath(parsed_args["T"], parsed_args["τ"],
        parsed_args["m"], parsed_args["m"],
        parsed_args["boltz_const"], parsed_args["sample_std"],
        parsed_args["D"], parsed_args["nu"], parsed_args["α"],
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
