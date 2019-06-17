using LatentOrderBookModel
using Test

include("../src/main.jl")

@testset "All Tests" begin

    @testset "Reproducibility" begin
        @test all(reaction_diffusion_path(45, 2300, 10, 238.745, 501, 2.0, 7.415,
            0.001, 5.0, 0.001, 1.0, 0.5) .== reaction_diffusion_path(45, 2300,
            10, 238.745, 501, 2.0, 7.415, 0.001, 5.0, 0.001, 1.0, 0.5))

        @test all(reaction_diffusion_path(51, 2300, 10, 238.745, 501, 2.0, 7.415,
            0.001, 5.0, 0.001, 1.0, 0.5) .== reaction_diffusion_path(51, 2300,
            10, 238.745, 501, 2, 7.415, 0.001, 5.0, 0.001, 1.0, 0.5))

        @test any(reaction_diffusion_path(45, 2300, 10, 238.745, 501, 2.0, 7.415,
            0.001, 5.0, 0.001, 1.0, 0.5) .!= reaction_diffusion_path(51, 2300,
            10, 238.745, 501, 2.0, 7.415, 0.001, 5.0, 0.001, 1.0, 0.5))
    end;

    @testset "Default Values" begin
        @test size(reaction_diffusion_path(),1) == 100
    end;

    @testset "Command Line Parse Default Arguments" begin
        @test parse_commandline() == Dict("μ" => 0.5,
            "T" => 100,
            "initial_mid_price" => 100.0,
            "λ" => 1.0,
            "σ" => 0.001,
            "SEED" => 1,
            "boltz_const" => 2.0,
            "η" => 0.001,
            "D" => 5.0,
            "n_spatial_points" => 101,
            "τ" => 10,
            "sample_std" => 4.0)
    end

    @testset "Main Default Arguments" begin
        @test all(main("return") .== reaction_diffusion_path())
    end
end
