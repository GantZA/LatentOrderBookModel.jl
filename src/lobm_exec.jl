module LOBMExec

using LatentOrderBookModel

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    parsed_args = parse_commandline()
    rdpp = ReactionDiffusionPricePath(parsed_args["T"], parsed_args["τ"],
        parsed_args["m"], parsed_args["m"],
        parsed_args["boltz_const"], parsed_args["sample_std"],
        parsed_args["D"], parsed_args["nu"], parsed_args["α"],
        SourceTerm(parsed_args["λ"], parsed_args["mu"]))
    print(rdpp(parsed_args["SEED"]))
    return 0
end


end
