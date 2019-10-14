function initial_conditions(rdpp::ReactionDiffusionPricePaths)
    Δx = (rdpp.x[end-1] - rdpp.x[2])/(rdpp.m)
    V₀ = rdpp.σ*rand(Normal(0.0, 1.0))

    # lower, middle, upper
    A = Tridiagonal(
        (V₀/(2.0*Δx) + rdpp.D/(Δx^2.0)) * ones(Float64, rdpp.m),
        ((-2.0*rdpp.D)/(Δx^2.0) - rdpp.nu) * ones(Float64, rdpp.m+1),
        (-V₀/(2.0*Δx) + rdpp.D/(Δx^2.0)) * ones(Float64, rdpp.m))

    A[1,1] = -rdpp.nu - rdpp.D/(Δx^2.0) + V₀/(2.0*Δx)
    A[end, end] = -rdpp.nu - rdpp.D/(Δx^2.0) - V₀/(2.0*Δx)

    B = .-[rdpp.source_term(xᵢ, rdpp.initial_mid_price) for xᵢ in rdpp.x[2:end-1]]

    ϕ = A \ B
    return ϕ
end

function extract_mid_price(lob_density)
    
end


function dtrw_solver(rdpp::ReactionDiffusionPricePaths)

    ϕ₀ = initial_conditions(rdpp)
    ϕ = ϕ₀[:]
    Φ = ones(Float64, rdpp.m+1, rdpp.T)
    Φ[:,1] = ϕ
    Δx = (rdpp.x[end] - rdpp.x[1])/(rdpp.m)
    Δt = (Δx^2) / (2.0*rdpp.D)
    # plot(1:501, ϕ)
    p =  ones(Float64, rdpp.T) * rdpp.initial_mid_price
    ϵ = rand(Normal(0.0,1.0), rdpp.T-1)
    P⁺s = ones(Float64, rdpp.m+1, rdpp.T-1)
    Ps = ones(Float64, rdpp.m+1, rdpp.T-1)
    P⁻s = ones(Float64, rdpp.m+1, rdpp.T-1)
    # Simulate SPDE
    @inbounds for n = 2:rdpp.T

        Vₜ = rdpp.σ*ϵ[n-1]
        V = -Vₜ.*rdpp.x./(2.0*rdpp.D)

        P⁺ = exp.(-rdpp.β.*V[3:end])
        P = exp.(-rdpp.β.*V[2:end-1])
        P⁻ = exp.(-rdpp.β.*V[1:end-2])
        Z = P⁺ .+ P .+ P⁻
        # Normalizing the probabilities
        P⁺ = P⁺ ./ Z
        P = P ./ Z
        P⁻ = P⁻ ./ Z


        P⁺s[:,n-1] = P⁺
        Ps[:,n-1] = P
        P⁻s[:,n-1] = P⁻

        ϕ₋₁ = ϕ₀[1]
        ϕₘ₊₁ = ϕ₀[end]

        ϕ[1] = P⁻[1] * ϕ₋₁ +
            P⁻[2] * ϕ₀[2] +
            P[1] * ϕ₀[1] -
            rdpp.nu * ϕ₀[1] +
            rdpp.source_term(rdpp.x[2], p[n-1])

        ϕ[end] = P⁺[end] * ϕₘ₊₁ +
            P⁺[end-1] * ϕ₀[end-1] +
            P[end] * ϕ₀[end] -
            rdpp.nu * ϕ₀[end] +
            rdpp.source_term(rdpp.x[end-1], p[n-1])

        # Compute Interior Points
        ϕ[2:end-1] = P⁺[1:end-2] .* ϕ₀[1:end-2] +
            P⁻[3:end] .* ϕ₀[3:end] +
            P[2:end-1] .* ϕ₀[2:end-1] -
            rdpp.nu * ϕ₀[2:end-1] +
            [rdpp.source_term(xᵢ, p[n-1]) for xᵢ in rdpp.x[3:end-2]]

        # The 'mass' at site 'j' at the next time step is the mass at 'j-1'
        # times the probability of right plus the mass at 'j+1' times the
        # probability of jumping left plus the mass at site 'j' that self
        # jumps back to site 'j'.

        ϕ₀ = ϕ[:]
        Φ[:,n] = ϕ[:]
        mid_price_ind = argmin(abs.(ϕ₀))
        p[n] = rdpp.x[mid_price_ind]
    end

    return Φ, p, P⁺s, Ps ,P⁻s
end
