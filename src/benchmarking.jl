using BenchmarkTools
using LatentOrderBookModel
@benchmark reaction_diffusion_path(45, 2300, 10, 238.745, 501, 2, 7.415,
    0.001, 5.0, 0.001, 1.0, 0.5)




function bench_rdp(seed)
    print(reaction_diffusion_path(seed, 2300, 10, 238.745, 501, 2, 7.415,
        0.001, 5.0, 0.001, 1.0, 0.5))
end

@benchmark bench_rdp.([565, 7362, 72382])
