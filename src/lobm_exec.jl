module LOBMExec

using LatentOrderBookModel

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    parsed_args = parse_commandline()
    rdpp = ReactionDiffusionPricePaths(
        1,
        parsed_args["T"],
        parsed_args["p₀"],
        parsed_args["M"],
        parsed_args["β"],
        parsed_args["L"],
        parsed_args["D"],
        parsed_args["σ"],
        parsed_args["nu"],
        SourceTerm(parsed_args["λ"], parsed_args["μ"]),
    )

    lob_densities,
    price_paths,
    mid_price_bars,
    P⁺s,
    P⁻s = rdpp(parsed_args["SEED"])
    print(mid_price_bars[:, 1])
    return 0
end


end
