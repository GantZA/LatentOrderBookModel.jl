using LatentOrderBookModel
using Test

include("../src/main.jl")

@testset "All Tests" begin
    rdpp = ReactionDiffusionPricePaths()
    @testset "Reproducibility" begin
        rdpp_1 = ReactionDiffusionPricePaths(1, 200, 238.745, 500, 2.0, 100,
            2.0, 0.01, 0.0, SourceTerm(1.0, 0.5))
        rdpp_2 = ReactionDiffusionPricePaths(1, 200, 238.745, 500, 2.0, 100,
            2.0, 0.01, 0.0, SourceTerm(1.0, 0.5))
        @test all(rdpp_1(45) .== rdpp_2(45))

        @test all(rdpp_1(51) .== rdpp_2(51))

        @test all(rdpp_1(4565756745) .== rdpp_2(4565756745))
    end;

    @testset "Default Values" begin
        lob_densities, price_paths, mid_price_bars, P⁺s, P⁻s = rdpp()
        @test size(mid_price_bars) == (101,1)
    end;

    @testset "Command Line Parse Default Arguments" begin
        @test parse_commandline() == Dict("μ" => 0.5,
            "num_paths" => 1,
            "T" => 100,
            "p₀" => 100.0,
            "λ" => 1.0,
            "SEED" => 1,
            "β" => 1.0,
            "nu" => 0.0,
            "σ" => 1.0,
            "D" => 4.0,
            "M" => 100,
            "L" => 50.0)
    end

    @testset "Main Default Arguments" begin
        @test all(main("return") .== rdpp(1))
    end
end
