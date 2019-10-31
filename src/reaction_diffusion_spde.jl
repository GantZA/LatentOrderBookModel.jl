function initial_conditions_steady_state(rdpp::ReactionDiffusionPricePaths)
    φ = [
            (xᵢ-rdpp.p₀) * rdpp.source_term.λ/(2 * rdpp.D * rdpp.source_term.μ) * exp((-rdpp.source_term.μ * rdpp.L^2) / 4) +
            (sqrt(pi) * rdpp.source_term.λ * erf(sqrt(rdpp.source_term.μ) * (rdpp.p₀ - xᵢ)))/(4 * rdpp.D * rdpp.source_term.μ^(3/2))
            for xᵢ in rdpp.x
        ]
    return φ
end


function extract_mid_price(rdpp, lob_density)
    mid_price_ind = 2
    while lob_density[mid_price_ind] > 0
        mid_price_ind += 1
    end
    x1, y1 = rdpp.x[mid_price_ind-1], lob_density[mid_price_ind-1]
    x2, y2 = rdpp.x[mid_price_ind], lob_density[mid_price_ind]
    m = (y1-y2)/(x1-x2)
    c = y1-m*x1
    mid_price = round(-c/m, digits=2)
    return mid_price
end


function dtrw_solver(rdpp::ReactionDiffusionPricePaths)

    φ = ones(Float64, rdpp.M+1, rdpp.T)
    φ[:,1] = initial_conditions_steady_state(rdpp)
    Δx = rdpp.L/rdpp.M
    Δt = (Δx^2) / (2.0*rdpp.D)

    p =  ones(Float64, rdpp.T) * rdpp.p₀
    ϵ = rand(Normal(0.0,1.0), rdpp.T-1)
    P⁺s = ones(Float64, rdpp.T-1)
    P⁻s = ones(Float64, rdpp.T-1)

    @inbounds for n = 1:rdpp.T-1
        Vₜ =  sign(ϵ[n]) * min(abs(rdpp.σ*ϵ[n]), Δx/Δt)
        V = (-Vₜ .* rdpp.x) ./ (2.0*rdpp.D)

        P⁺ = 1/(1 + exp(-rdpp.β * Vₜ * Δx / rdpp.D))
        P⁻ = 1 - P⁺

        P⁺s[n] .= P⁺
        P⁻s[n] .= P⁻

        φ₋₁ = φ[1,n] * (1-(Vₜ*Δx)/(2*rdpp.D))
        φₘ₊₁ = φ[end,n] * (1+(Vₜ*Δx)/(2*rdpp.D))

        φ[1,n+1] = P⁺ * φ₋₁ +
            P⁻ * φ[2,n] +
            rdpp.source_term(rdpp.x[1], p[n])

        φ[end,n+1] = P⁻ * φₘ₊₁ +
            P⁺ * φ[end-1,n] +
            rdpp.source_term(rdpp.x[end], p[n])

        φ[2:end-1,n+1] = P⁺ * φ[1:end-2,n] +
            P⁻ * φ[3:end,n] +
            [rdpp.source_term(xᵢ, p[n]) for xᵢ in rdpp.x[2:end-1]]

        p[n+1] = extract_mid_price(rdpp, φ[:,n+1])
    end

    return φ, p, P⁺s ,P⁻s
end
