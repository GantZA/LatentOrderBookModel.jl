using LatentOrderBookModel
using Test

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

end
