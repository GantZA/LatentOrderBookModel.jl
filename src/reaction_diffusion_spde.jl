function initial_conditions_numerical(rdpp::ReactionDiffusionPricePaths, pₙ)
    ϵ = rand(Normal(0.0, 1.0))
    V₀ = sign(ϵ) * min(abs(rdpp.σ * ϵ), rdpp.Δx / rdpp.Δt)

    A = Tridiagonal(
        (V₀/(2.0*rdpp.Δx) + rdpp.D/(rdpp.Δx^2)) * ones(Float64, rdpp.M),
        ((-2.0*rdpp.D)/(rdpp.Δx^2) - rdpp.nu) * ones(Float64, rdpp.M+1),
        (-V₀/(2.0*rdpp.Δx) + rdpp.D/(rdpp.Δx^2)) * ones(Float64, rdpp.M))

    A[1,1] = (-rdpp.D)/(rdpp.Δx^2) - rdpp.nu - V₀/(2.0*rdpp.Δx)
    A[end, end] = (-rdpp.D)/(rdpp.Δx^2) - rdpp.nu - V₀/(2.0*rdpp.Δx)

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
    while (lob_density[mid_price_ind] > 0) | (lob_density[mid_price_ind+1]>lob_density[mid_price_ind])
        mid_price_ind += 1
    end

    y1 = lob_density[mid_price_ind-1]
    y2 = lob_density[mid_price_ind]
    x1 = rdpp.x[mid_price_ind-1]

    mid_price = round(-(y1 * rdpp.Δx)/(y2 - y1) + x1, digits = 2)
    return mid_price
end


function calculate_right_jump_probability(Z)
    if Z > 10
        return 1.0
    elseif Z < -10
        return 0.0
    else
        return (exp(Z))/(exp(-Z) + exp(Z) + 1)
    end
end

function calculate_left_jump_probability(Z)
    if Z > 10
        return 0.0
    elseif Z < -10
        return 1.0
    else
        return (exp(-Z))/(exp(-Z) + exp(Z) + 1)
    end
end

function calculate_self_jump_probability(Z)
    if Z > 10
        return 0.0
    elseif Z < -10
        return 0.0
    else
        return 1/(exp(-Z) + exp(Z) + 1)
    end
end

function calculate_jump_probabilities(rdpp, Vₜ)
    Z = (rdpp.β * Vₜ * rdpp.Δx) / (2* rdpp.D)
    return calculate_right_jump_probability(Z), calculate_left_jump_probability(Z), calculate_self_jump_probability(Z)
end

function get_sub_period_time(rdpp, t, time_steps)
    τ = rand(Exponential(1/rdpp.α))
    remaining_time = time_steps - t
    τ_periods = min(floor(Int, τ/rdpp.Δt), remaining_time)
    return τ, τ_periods
end


function intra_time_period_simulate(rdpp, φ, p)
    ϵ = rand(Normal(0.0, 1.0))
    Vₜ = sign(ϵ) * min(abs(rdpp.σ * ϵ), rdpp.Δx / rdpp.Δt)

    P⁺, P⁻, P = calculate_jump_probabilities(rdpp, Vₜ)

    φ₋₁ = φ[1]
    φₘ₊₁ = φ[end]
    φ_next = zeros(Float64, size(φ,1))

    φ_next[1] = P⁺ * φ₋₁ + P⁻ * φ[2] + rdpp.source_term(rdpp.x[1], p)
    φ_next[end] = P⁻ * φₘ₊₁ + P⁺ * φ[end-1] + rdpp.source_term(rdpp.x[end], p)
    φ_next[2:end-1] = P⁺ * φ[1:end-2] + P⁻ * φ[3:end] +
        [rdpp.source_term(xᵢ, p) for xᵢ in rdpp.x[2:end-1]]

    return φ_next, P⁺, P⁻, P
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
    Ps = ones(Float64, time_steps)

    t = 1
    calendar_time_index = 1
    while t < time_steps
        τ, τ_periods = get_sub_period_time(rdpp, t, time_steps)
        φ[:, t] = initial_conditions_numerical(rdpp, p[t])
        for τₖ = 1:τ_periods
            φ[:, t+1], P⁺s[t], P⁻s[t], Ps[t]  = intra_time_period_simulate(rdpp, φ[:, t], p[calendar_time_index])
            p[t+1] = extract_mid_price(rdpp, φ[:, t+1])
            if t*rdpp.Δt > calendar_time_index
                calendar_time_index += 1
                mid_prices[calendar_time_index] = extract_mid_price(rdpp, φ[:, t+1])
            end
            t += 1
        end
    end

    for t = 1:
        τ = floor(Int, (1 + mod(n-1, rdpp.Δt))/rdpp.Δt)
        φ[:, initial_lob_index] = initial_conditions_numerical(rdpp, mid_prices[n])

        for t = 1:τ
            running_index = initial_lob_index + t
            ϵ = rand(Normal(0.0, 1.0))
            Vₜ = sign(ϵ) * min(abs(rdpp.σ * ϵ), rdpp.Δx / rdpp.Δt)

            P⁺, P⁻, P = calculate_jump_probabilities(rdpp, Vₜ)

            P⁺s[running_index-1] = P⁺
            P⁻s[running_index-1] = P⁻
            P[running_index-1] = P

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
        initial_lob_index += τ + 1
        mid_prices[n+1] = p[initial_lob_index-1]
    end
    return φ, p, mid_prices, P⁺s, P⁻s
end
