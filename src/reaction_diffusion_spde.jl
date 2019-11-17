function initial_conditions_numerical(rdpp::ReactionDiffusionPricePaths, pₙ)
    ϵ = rand(Normal(0.0, 1.0))
    V₀ = sign(ϵ) * min(abs(rdpp.σ * ϵ), rdpp.Δx / rdpp.Δt)

    A = Tridiagonal(
        (V₀/(2.0*rdpp.Δx) + rdpp.D/(rdpp.Δx^2)) * ones(Float64, rdpp.M),
        ((-2.0*rdpp.D)/(rdpp.Δx^2) - rdpp.nu) * ones(Float64, rdpp.M+1),
        (-V₀/(2.0*rdpp.Δx) + rdpp.D/(rdpp.Δx^2)) * ones(Float64, rdpp.M))

    A[1,1] = (-rdpp.D)/(rdpp.Δx^2) - rdpp.nu + V₀/(2.0*rdpp.Δx)
    A[end, end] = 2.0*rdpp.D/rdpp.Δx^2

    B = .-[rdpp.source_term(xᵢ, pₙ) for xᵢ in rdpp.x]
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
    mid_prices = ones(Float64, rdpp.T + 1)

    p[1] = rdpp.p₀
    mid_prices[1] = rdpp.p₀

    P⁺s = ones(Float64, time_steps)
    P⁻s = ones(Float64, time_steps)

    initial_lob_index = 1

    @inbounds for n = 1:rdpp.T
        τ = floor(Int, (1 + mod(n-1, rdpp.Δt))/rdpp.Δt)
        φ[:, initial_lob_index] = initial_conditions_numerical(rdpp, mid_prices[n-1])

        for t = 1:τ
            running_index = initial_lob_index + t
            ϵ = rand(Normal(0.0, 1.0))
            Vₜ = sign(ϵ) * min(abs(rdpp.σ * ϵ), rdpp.Δx / rdpp.Δt)

            P⁺ = 1 / (1 + exp(-rdpp.β * Vₜ * rdpp.Δx / rdpp.D))
            P⁻ = 1 - P⁺

            P⁺s[running_index-1] = P⁺
            P⁻s[running_index-1] = P⁻

            φ₋₁ = φ[1, running_index-1]
            φₘ₊₁ = φ[end, running_index-1]

            φ[1, running_index] = P⁺ * φ₋₁ + P⁻ * φ[2, running_index-1] +
                rdpp.source_term(rdpp.x[1], p[running_index-1])

            φ[end, running_index] = P⁻ * φₘ₊₁ + P⁺ * φ[end-1, running_index-1] +
                          rdpp.source_term(rdpp.x[end], p[running_index-1])

            φ[2:end-1, running_index] = P⁺ * φ[1:end-2, running_index-1] +
                P⁻ * φ[3:end, running_index-1] +
                [rdpp.source_term(xᵢ, p[running_index-1]) for xᵢ in rdpp.x[2:end-1]]

            p[running_index] = extract_mid_price(rdpp, φ[:, running_index])
        end
        initial_lob_index += initial_lob_index + τ + 1
        mid_prices[n+1] = p[initial_lob_index-1]
    end
    return φ, p, mid_prices, P⁺s, P⁻s
end
