using Statistics
using Distributions
using Random

# Generate a time series of stock prices by repeatedly solving the PDE
# from Mastromatteo et al. (2014) to model the latent order book.
# using the previously generated mid price as the new initial mid price and
# boundary mid point during each iteration to generate a sequence of T prices,
# where T is the number of measured prices stored in the loaded
# ObjectiveFunction object.

function (rdpp::ReactionDiffusionPricePath)(seed::Int=-1)
    if seed == -1
            Random.seed!()
        else
            Random.seed!(seed)
    end
    p = rdpp.initial_mid_price * ones(Float64, rdpp.T)
    for t in 1:rdpp.T-1
        a = max(0.0, p[t] - 3.0*rdpp.sample_std)
        b = p[t] + 3.0*rdpp.sample_std

        # Price Grid
        x = collect(Float64, range(a, b, length=rdpp.n_spatial_points))

        dtrw_solver = DTRWSolver(a, b, x, p[t])
        U = dtrw_solver(rdpp)

        mid_price_ind = argmin(abs.(U))
        p[t+1] = x[mid_price_ind]
    end
    return p
end
