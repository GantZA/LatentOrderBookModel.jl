using LatentOrderBookModel
using Test

include("../src/main.jl")

@testset "All Tests" begin
    rdpp = ReactionDiffusionPricePaths()
    @testset "Reproducibility" begin
        rdpp_1 = ReactionDiffusionPricePaths(1, 1000, 238.745, 500, 2.0, 7.415,
            5.0, 0.001, 0.0, SourceTerm(1.0, 0.5))
        rdpp_2 = ReactionDiffusionPricePaths(1, 1000, 238.745, 500, 2.0, 7.415,
            5.0, 0.001, 0.0, SourceTerm(1.0, 0.5))
        @test all(rdpp_1(45) .== rdpp_2(45))

        @test all(rdpp_1(51) .== rdpp_2(51))

        @test all(rdpp_1(4565756745) .== rdpp_2(4565756745))
    end;

    @testset "Default Values" begin
        lob_densities, price_paths, P⁺s, P⁻s = rdpp()
        @test size(price_paths) == (100,1)
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
            "σ" => 0.001,
            "D" => 5.0,
            "M" => 100,
            "sample_std" => 4.0)
    end

    @testset "Main Default Arguments" begin
        @test all(main("return") .== rdpp(1))
    end
end
