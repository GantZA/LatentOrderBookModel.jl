using LatentOrderBookModel
using Test

include("../src/main.jl")

@testset "All Tests" begin
    rdpp = ReactionDiffusionPricePath()
    @testset "Reproducibility" begin
        rdpp_1 = ReactionDiffusionPricePath(2300, 10, 238.745, 501, 2.0, 7.415,
            5.0, 0.001, 1.0, SourceTerm(1.0, 0.5))
        rdpp_2 = ReactionDiffusionPricePath(2300, 10, 238.745, 501, 2.0, 7.415,
            5.0, 0.001, 1.0, SourceTerm(1.0, 0.5))
        @test all(rdpp_1(45) .== rdpp_2(45))

        @test all(rdpp_1(51) .== rdpp_2(51))

        @test all(rdpp_1(4565756745) .== rdpp_2(4565756745))
    end;

    @testset "Default Values" begin
        @test size(rdpp(),1) == 100
    end;

    @testset "Command Line Parse Default Arguments" begin
        @test parse_commandline() == Dict("μ" => 0.5,
            "T" => 100,
            "m" => 100.0,
            "λ" => 1.0,
            "SEED" => 1,
            "boltz_const" => 1.0,
            "nu" => 0.001,
            "α" => 1.0,
            "D" => 5.0,
            "m" => 101,
            "τ" => 10,
            "sample_std" => 4.0)
    end

    @testset "Main Default Arguments" begin

        @test all(main("return") .== rdpp(1))
    end
end
