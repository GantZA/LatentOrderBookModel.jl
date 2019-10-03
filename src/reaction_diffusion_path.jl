
mutable struct ReactionDiffusionPricePath
    T::Int64
    initial_mid_price::Float64
    m::Int64
    β::Float64
    sample_std::Float64
    D::Float64
    nu::Float64
    σ::Float64
    α::Float64
    source_term::SourceTerm
end

function (rdpp::ReactionDiffusionPricePath)(seed::Int=-1)
    if seed == -1
            Random.seed!()
        else
            Random.seed!(seed)
    end
    p = rdpp.initial_mid_price * ones(Float64, rdpp.T)
    # lob_densities = zeros(Float64, rdpp.m+1, rdpp.T)
    for k in 1:rdpp.T-1
        p[k+1] = dtrw_solver(rdpp, p[k])
    end
    return p
end

ReactionDiffusionPricePath(dict)=ReactionDiffusionPricePath(
    dict["T"], dict["initial_mid_price"], dict["m"],
    dict["β"], dict["sample_std"], dict["D"], dict["nu"],
    dict["σ"],dict["α"], SourceTerm(dict["λ"], dict["μ"]))

ReactionDiffusionPricePath(;
    T::Int64=100,
    initial_mid_price::Real=100.0,
    m::Int64=100,
    β::Real=1.0,
    sample_std::Real=4.0,
    D::Real=5.0,
    nu::Real=0.001,
    σ::Real=0.001,
    α::Real=1.0,
    λ::Real=1.0,
    μ::Real=0.5) =
    ReactionDiffusionPricePath(T, initial_mid_price, m, β, sample_std,
        D, nu, σ, α, SourceTerm(λ, μ))
