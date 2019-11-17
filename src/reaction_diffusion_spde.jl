function initial_conditions_numerical(rdpp::ReactionDiffusionPricePaths, pₙ)

    ϵ = rand(Normal(0.0, 1.0))
    V₀ = sign(ϵ) * min(abs(rdpp.σ * ϵ), rdpp.Δx / rdpp.Δt)

    A = Tridiagonal(
        (V₀/(2.0*rdpp.Δx) + rdpp.D/(rdpp.Δx^2.0)) * ones(Float64, rdpp.M),
        ((-2.0*rdpp.D)/(rdpp.Δx^2.0) - rdpp.nu) * ones(Float64, rdpp.M+1),
        (-V₀/(2.0*rdpp.Δx) + rdpp.D/(rdpp.Δx^2.0)) * ones(Float64, rdpp.M))

    A[1,2] = 2.0*rdpp.D/rdpp.Δx^2
    A[end, end-1] = 2.0*rdpp.D/rdpp.Δx^2

    B = .-[rdpp.source_term(xᵢ, pₙ) for xᵢ in rdpp.x]= A \ B
    φ = A \ B
    return φ
end


function sample_mid_price_path(rdpp, Δt, price_path)
    mid_prices = zeros(Float64, rdpp.T + 1)
    mid_prices[1] = price_path[1]
    for t = 1:rdpp.T
        close_ind = floor(Int, t / Δt)
        mid_prices[t+1] = price_path[close_ind]
    end
    return mid_prices
end


function extract_mid_price(rdpp, lob_density)
    mid_price_ind = 2
    while lob_density[mid_price_ind] > 0
        mid_price_ind += 1
    end

    y1 = lob_density[mid_price_ind-1]
    y2 = lob_density[mid_price_ind]
    x1 = rdpp.x[mid_price_ind-1]

    mid_price = round(-(y1 * rdpp.Δx)/(y2 - y1) + x1, digits = 2)
    return mid_price
end


function dtrw_solver(rdpp::ReactionDiffusionPricePaths)
    time_steps = floor(Int, rdpp.T / rdpp.Δt)

    φ = ones(Float64, rdpp.M + 1, time_steps + 1)



    p = ones(Float64, time_steps + 1)
    mid_prices = ones(Float64, rdpp.T)

    p[1] = rdpp.p₀
    mid_prices[1] = rdpp.p₀

    ϵ = rand(Normal(0.0, 1.0), time_steps)
    P⁺s = ones(Float64, time_steps)
    P⁻s = ones(Float64, time_steps)

    @inbounds for n = 2:rdpp.T
        τ = floor(Int, (1 + mod(n-1, rdpp.Δt))/rdpp.Δt)
        φ[:, n-1] = initial_conditions_numerical(rdpp, mid_prices[n-1])
        for t = 1:τ

        end



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
