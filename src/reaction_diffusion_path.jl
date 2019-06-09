using Statistics
using Distributions
using Random

include("reaction_diffusion_spde.jl")

# Generate a time series of stock prices by repeatedly solving the PDE
# from Mastromatteo et al. (2014) to model the latent order book.
# using the previously generated mid price as the new initial mid price and
# boundary mid point during each iteration to generate a sequence of T prices,
# where T is the number of measured prices stored in the loaded
# ObjectiveFunction object.


mutable struct RDP_Params
    SEED::Int64
    T::Int64
    τ::Int64
    initial_mid_price::Float64
    n_spatial_points::Int64
    boltz_const::Float64
    sample_std::Float64
    σ::Float64
    D::Float64
    η::Float64
end

function reaction_diffusion_path(rdp_params, st_params)

    Random.seed!(rdp_params.SEED)

    # get reaction diffusion path parameters from rdp_params
    p0          = rdp_params.initial_mid_price
    T           = rdp_params.T
    sample_std  = rdp_params.sample_std
    σ           = rdp_params.σ
    β           = rdp_params.boltz_const
    n           = rdp_params.n_spatial_points
    D           = rdp_params.D
    η           = rdp_params.η
    τ           = rdp_params.τ

    # construct vector of simulated prices, filling with p0
    p = p0 * ones(T)

    for t in 1:T-1
        # Boundaries of prices [a,b]
        a = max(0, p[t] - 3*sample_std)
        b = p[t] + 3*sample_std

        # Price Grid
        x = collect(range(a, b, length=n))

        # Generate the Advective Coefficient using a seed for reproducibility
        V₀ = rand(Normal(0, σ))

        # Create struct for DTRW extra params
        dtrw_rd_params = DTRW_RD_Params([a, b], x, V₀, p[t])

        # Iterate DTRW Solver τ time steps
        U =  DTRW_reaction_diffusion(dtrw_rd_params, rdp_params,
            st_params)

        # Extract Mid-Price by finding the index where the density functions
        # meet. The price at that index is the next simulated midprice
        mid_price_ind = argmin(abs.(U))
        p[t+1] = x[mid_price_ind]
    end
    return p
end
