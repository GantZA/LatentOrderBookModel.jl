function initial_conditions(rdpp::ReactionDiffusionPricePaths)
    Δx = (rdpp.x[end] - rdpp.x[1])/(rdpp.m)


    # lower, middle, upper
    ϕ = [rdpp.source_term.λ * tanh(rdpp.source_term.μ *
        (rdpp.initial_mid_price -(i+1/2)*Δx)) for i in 1:rdpp.m+1]
    return ϕ
end


function dtrw_solver(rdpp::ReactionDiffusionPricePaths)

    ϕ₀ = initial_conditions(rdpp)
    ϕ = ϕ₀[:]
    Φ = ones(Float64, rdpp.m+1, rdpp.T)
    Φ[:,1] = ϕ
    Δx = (rdpp.x[end] - rdpp.x[1])/(rdpp.m)
    Δt = rdpp.r * (Δx^2) / (2.0*rdpp.D)
    # plot(1:501, ϕ)
    p =  ones(Float64, rdpp.T) * rdpp.initial_mid_price
    ϵ = rand(Normal(0.0,2.0), rdpp.T-1)
    # P⁺s = ones(Float64, rdpp.m+1, rdpp.T-1)
    # Ps = ones(Float64, rdpp.m+1, rdpp.T-1)
    # P⁻s = ones(Float64, rdpp.m+1, rdpp.T-1)
    # Simulate SPDE
    @inbounds for n = 2:rdpp.T


        Vₜ = ϵ[n-1]

        P⁺ = (exp(-rdpp.nu * Δt))/(1+exp(-Vₜ/rdpp.D * Δx))*rdpp.r
        P⁻ = (exp(-rdpp.nu * Δt))/(1+exp(Vₜ/rdpp.D * Δx))*rdpp.r
        ϕ[2:end-1] = P⁺ * ϕ₀[1:end-2] +
            P⁻ * ϕ₀[3:end] +
            (1-rdpp.r) * exp(-rdpp.nu*Δt) * ϕ₀[2:end-1] -
            (1-exp(-rdpp.nu*Δt)) * ϕ₀[2:end-1] +
            [rdpp.source_term(i, Δx, p[n-1]) for i in 2:rdpp.m] * Δt


        ϕ[1] = P⁺ * ϕ₀[1] +
            P⁻ * ϕ₀[2] +
            (1-rdpp.r) * exp(-rdpp.nu*Δt) * ϕ₀[1] -
            (1-exp(-rdpp.nu*Δt)) * ϕ₀[1] +
            rdpp.source_term(1, Δx, p[n-1]) * Δt

        ϕ[end] = P⁺ * ϕ₀[end-1] +
            P⁻ * ϕ₀[end] +
            (1-rdpp.r) * exp(-rdpp.nu*Δt) * ϕ₀[end] -
            (1-exp(-rdpp.nu*Δt)) * ϕ₀[end] +
            rdpp.source_term(rdpp.m, Δx, p[n-1]) * Δt

        ϕ₀ = ϕ[:]
        Φ[:,n] = ϕ[:]
        mid_price_ind = argmin(abs.(ϕ₀))

        p[n] =  (ϕ[mid_price_ind] * Δx)/(ϕ[mid_price_ind]-ϕ[mid_price_ind+1]) +
            (mid_price_ind+1/2)*Δx
    end

    return p, Φ
end
