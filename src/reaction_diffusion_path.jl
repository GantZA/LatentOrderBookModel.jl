# Generate a time series of stock prices by repeatedly solving the SPDE
# from Mastromatteo et al. (2014) to model the latent order book.
# using the previously generated mid price as the new initial mid price and
# boundary mid point during each iteration to generate a sequence of T prices,
# where T is the number of measured prices stored in the loaded
# ObjectiveFunction object.

mutable struct ReactionDiffusionPricePaths
    num_paths::Int64
    T::Int64
    p₀::Float64
    M::Int64
    β::Float64
    sample_std::Float64
    D::Float64
    σ::Float64
    source_term::SourceTerm
    x::Array{Float64, 1}
end
function ReactionDiffusionPricePaths(num_paths, T, p₀, M, β,
    sample_std, D, σ, source_term)
    x₀ = max(0.0, p₀ - 3.0*sample_std)
    xₘ = p₀ + 3.0*sample_std
    x = collect(Float64, range(x₀, xₘ, length=M+1))
    return ReactionDiffusionPricePaths(num_paths, T, p₀, M, β,
        sample_std, D, σ, source_term, x)
end


function (rdpp::ReactionDiffusionPricePaths)(seed::Int=-1)
    if seed == -1
            seeds = Int.(rand(MersenneTwister(), UInt32, rdpp.num_paths))
        else
            seeds = Int.(rand(MersenneTwister(seed), UInt32, rdpp.num_paths))
    end

    price_paths = ones(Float64, rdpp.T, rdpp.num_paths) * rdpp.p₀
    lob_densities = zeros(Float64, rdpp.M+1, rdpp.T, rdpp.num_paths)
    P⁺s = ones(Float64, rdpp.M+1, rdpp.T-1, rdpp.num_paths)
    # Ps = ones(Float64, rdpp.M+1, rdpp.T-1, rdpp.num_paths)
    P⁻s = ones(Float64, rdpp.M+1, rdpp.T-1, rdpp.num_paths)

    for path in 1:rdpp.num_paths
         Random.seed!(seeds[path])
        lob_densities[:, :, path], price_paths[:, path], P⁺s[:,:,path],
            P⁻s[:,:,path]  = dtrw_solver(rdpp)
    end

    return lob_densities, price_paths, P⁺s, P⁻s
end


ReactionDiffusionPricePaths(dict)=ReactionDiffusionPricePaths(
    dict["num_paths"], dict["T"], dict["p₀"], dict["M"],
    dict["boltz_const"], dict["sample_std"], dict["D"],
    dict["σ"], SourceTerm(dict["λ"], dict["μ"]))


ReactionDiffusionPricePaths(;num_paths=1,T::Int64=100,
    p₀::Real=100.0, M::Int64=100, boltz_const::Real=1.0,
    sample_std::Real=4.0, D::Real=5.0, σ::Real=1.0,
    λ::Real=1.0, μ::Real=0.5) =
    ReactionDiffusionPricePaths(num_paths, T, p₀,
    M, boltz_const, sample_std, D, σ,
    SourceTerm(λ, μ))
