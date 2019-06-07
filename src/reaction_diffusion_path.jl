using Statistics
using Distributions

include("reaction_diffusion_spde.jl")
include("mom_basic.jl")

# Generate a time series of stock prices by repeatedly solving the PDE
# from Mastromatteo et al. (2014) to model the latent order book.
# using the previously generated mid price as the new initial mid price and
# boundary mid point during each iteration to generate a sequence of T prices,
# where T is the number of measured prices stored in the loaded
# ObjectiveFunction object.



function reaction_diffusion_path(bar_data, σ, D, ν, source_func::Function,
    source_params)

    T = size(bar_data,1)
    # determine initial mid-price
    prices = bar_data.price

    # p0 = mean(prices)
    p0 = prices[1]

    # construct price vector
    p = p0 * ones(T)

    # Number of Spatial (price) discretization points
    n = 501

    # Boltzmann Constant
    β = 2

    sample_std = std(prices)

    for t in 1:T-1
        # Boundaries of prices [a,b]
        a = max(0,p[t] - 3*sample_std)
        b = p[t] + 3*sample_std

        # Price Grid
        x = collect(range(a,b, length= n))

        V₀ = rand(Normal(0,σ))

        # Iterate Solver
        U =  DTRW_reaction_diffusion(a,b,D,V₀,ν,source_params,β,p[t],n,10)

        # Extract Mid-Price
        mid_price_ind = argmin(abs.(U))
        p[t+1] = x[mid_price_ind]
    end
    return p
end

# Params: log_prices, σ, D, ν, source_func::Function,source_params
