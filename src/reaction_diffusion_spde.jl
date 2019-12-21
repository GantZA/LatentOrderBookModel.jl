function initial_conditions_steady_state(rdpp::ReactionDiffusionPricePaths)
    φ = [(xᵢ - rdpp.p₀) * rdpp.source_term.λ /
         (2 * rdpp.D * rdpp.source_term.μ) *
         exp((-rdpp.source_term.μ * rdpp.L^2) / 4) +
         (sqrt(pi) * rdpp.source_term.λ *
          erf(sqrt(rdpp.source_term.μ) * (rdpp.p₀ - xᵢ))) /
         (4 * rdpp.D * rdpp.source_term.μ^(3 / 2)) for xᵢ in rdpp.x]
    return φ
end


function sample_mid_price_path(rdpp, Δt, price_path)
    mid_prices = zeros(Float64, rdpp.T + 1)
    mid_prices[1] = price_path[1]
    for t = 1:rdpp.T
        close_ind = floor(Int, t / Δt) + 1
        mid_prices[t+1] = price_path[close_ind]
    end
    return mid_prices
end


function extract_mid_price(rdpp, lob_density)
    mid_price_ind = 2
    while (lob_density[mid_price_ind] > 0) | (lob_density[mid_price_ind+1]>lob_density[mid_price_ind])
        mid_price_ind += 1
        if mid_price_ind > rdpp.M
            mid_price_ind += 1
            break
        end
    end
    y1 = lob_density[mid_price_ind-1]
    y2 = lob_density[mid_price_ind]
    x1 = rdpp.x[mid_price_ind-1]
    mid_price = round(-(y1 * rdpp.Δx)/(y2 - y1) + x1, digits = 2)
    return mid_price
end


function dtrw_solver(rdpp::ReactionDiffusionPricePaths)
    Δx = rdpp.L / rdpp.M
    Δt = (Δx^2) / (2.0 * rdpp.D)
    time_steps = floor(Int, rdpp.T / rdpp.Δt)

    φ = ones(Float64, rdpp.M + 1, time_steps + 1)
    φ[:, 1] = initial_conditions_steady_state(rdpp)


    p = ones(Float64, time_steps + 1)
    p[1] = rdpp.p₀
    ϵ = rand(Normal(0.0, 1.0), time_steps)
    P⁺s = ones(Float64, time_steps)
    P⁻s = ones(Float64, time_steps)

    @inbounds for n = 1:time_steps
        Vₜ = sign(ϵ[n]) * min(abs(rdpp.σ * ϵ[n]), rdpp.Δx / rdpp.Δt)
        V = (-Vₜ .* rdpp.x) ./ (2.0 * rdpp.D)

        P⁺ = 1 / (1 + exp(-rdpp.β * Vₜ * rdpp.Δx / rdpp.D))
        P⁻ = 1 - P⁺

        P⁺s[n] = P⁺
        P⁻s[n] = P⁻

        φ₋₁ = (1-Vₜ*rdpp.Δx/rdpp.D) * φ[1, n]
        φₘ₊₁ = (1+Vₜ*rdpp.Δx/rdpp.D) * φ[end, n]

        φ[1, n+1] = P⁺ * φ₋₁ + P⁻ * φ[2, n] + rdpp.source_term(rdpp.x[1], p[n])

        φ[end, n+1] = P⁻ * φₘ₊₁ + P⁺ * φ[end-1, n] +
                      rdpp.source_term(rdpp.x[end], p[n])

        φ[2:end-1, n+1] = P⁺ * φ[1:end-2, n] + P⁻ * φ[3:end, n] +
                          [rdpp.source_term(xᵢ, p[n]) for xᵢ in rdpp.x[2:end-1]]

        p[n+1] = extract_mid_price(rdpp, φ[:, n+1])
    end
    mid_price_bars_close = sample_mid_price_path(rdpp, rdpp.Δt, p)
    return φ, p, mid_price_bars_close, P⁺s, P⁻s
end
