
function initial_conditions(rdpp::ReactionDiffusionPricePath, x, p_n)
    Δx = (x[end] - x[1])/rdpp.m
    V₀ = rand(Normal(0.0, rdpp.σ))

    A = Tridiagonal(
        (V₀/(2.0*Δx) + rdpp.D/(Δx^2.0)) * ones(Float64, rdpp.m),
        ((-2.0*rdpp.D)/(Δx^2.0) - rdpp.nu) * ones(Float64, rdpp.m+1),
        (-V₀/(2.0*Δx) + rdpp.D/(Δx^2.0)) * ones(Float64, rdpp.m))

    A[1,2] = 2.0*rdpp.D/Δx^2
    A[end, end-1] = 2.0*rdpp.D/Δx^2

    B = .-[rdpp.source_term(xᵢ, p_n) for xᵢ in x]

    ϕ = A \ B
    return ϕ
end

function dtrw_solver(rdpp::ReactionDiffusionPricePath, p_n)

    p_new = p_n
    x₀ = max(0.0, p_n - 3.0*rdpp.sample_std)
    xₘ = p_n + 3.0*rdpp.sample_std
    x = collect(Float64, range(x₀, xₘ, length=rdpp.m+1))

    ϕ₀ = initial_conditions(rdpp, x, p_n)
    ϕ = ϕ₀[:]

    Δx = (x[end] - x[1])/(rdpp.m)
    Δt = (Δx^2) / (2.0*rdpp.D)

    τ = rand(Exponential(1/rdpp.α))
    N = max(1, ceil(Int64, τ/Δt))
    ϵ = rand(Normal(0.0,rdpp.σ), N)

    @inbounds for n = 1:N
        Vₜ = ϵ[n]
        V = -Vₜ*x/(2.0*rdpp.D)

        P⁺ = vcat(exp.(-rdpp.β*V[2:end]),exp(-rdpp.β*V[end]))
        P = exp.(-rdpp.β*V)
        P⁻ = vcat(exp(-rdpp.β*V[1]),exp.(-rdpp.β*V[1:end-1]))
        Z = P⁺ .+ P .+ P⁻

        P⁺ = P⁺ ./ Z
        P = P ./ Z
        P⁻ = P⁻ ./ Z

        ϕ[1] = P⁻[1] * ϕ₀[1] +
            P⁻[2] * ϕ₀[2] +
            P[1] * ϕ₀[1] -
            rdpp.nu * ϕ₀[1] +
            rdpp.source_term(x[1], p_new)

        ϕ[2:end-1] = P⁺[1:end-2] .* ϕ₀[1:end-2] .+
            P⁻[3:end] .* ϕ₀[3:end] .+
            P[2:end-1] .* ϕ₀[2:end-1] .-
            rdpp.nu * ϕ₀[2:end-1] .+
            [rdpp.source_term(xᵢ, p_new) for xᵢ in x[2:end-1]]

        ϕ[end] = P⁺[end-1] * ϕ₀[end-1] +
            P⁺[end] * ϕ₀[end] +
            P[end] * ϕ₀[end] -
            rdpp.nu * ϕ₀[end] +
            rdpp.source_term(x[end], p_new)

        ϕ₀ = ϕ[:]
        mid_price_ind = argmin(abs.(ϕ₀))
        p_new = x[mid_price_ind]
    end
    return p_new
end
