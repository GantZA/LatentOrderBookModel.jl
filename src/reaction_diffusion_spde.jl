using LinearAlgebra

# SPDE
# ϱₜ = -Vₜ * ϱₓ + D * ϱₓₓ - ν * ϱ(x,t) + λ * tanh(μ * (p(t) - x))

function source_term(x, source_params, p)
    λ = source_params[1]
    μ = source_params[2]
    return λ*tanh(μ*(p-x))
end


function stochastic_drift(α₀, α, ε)
    return α₀ + α*ε
end

function get_initial_conditions(a,b,Vₜ,D, ν, source_params, p, n)
    xs = collect(range(a, b, length = n))
    Δx = (b - a)/(n-1)
    A = Tridiagonal(
        (Vₜ/(2*Δx) + D/(Δx^2)) * ones(n-1),
        ((-2*D)/(Δx^2) - ν) * ones(n),
        (-Vₜ/(2*Δx) + D/(Δx^2)) * ones(n-1))

    A[1,2] = 2*D/Δx^2
    A[end, end-1] = 2*D/Δx^2

    B = .-[source_term(x, source_params,p) for x in xs]

    U = A \ B
    return U
end

function DTRW_reaction_diffusion(a,b,D,V₀, ν, source_params, β, p0, n, T)
    #ϱₜ = -Vₜ * ϱₓ + D * ϱₓₓ - ν * ϱ(x,t) + λ * tanh(μ * (p(t) - x))

    u0 = get_initial_conditions(a,b,V₀,D, ν, source_params, p0, n)
    # plt = plot(1:200, u0, label = "u0");
    u = u0[:]
    x = collect(range(a, b, length = n))
    Δx = (b - a)/(n - 1)
    Δt = (Δx^2) / (2*D)
    V = -V₀*x/(2*D*β)
    # Simulate PDE
    for i = 1:T
        # ϵ = rand(Normal(0,1))
        # vₜ = stochastic_drift(0.0,1.0,ϵ)
        jump_prob_right = vcat(exp.(-β*V[2:end]),exp(-β*V[end]))
        jump_prob_self = exp.(-β*V)
        jump_prob_left = vcat(exp(-β*V[1]),exp.(-β*V[1:end-1]))
        jump_prob_denominator = jump_prob_right + jump_prob_self + jump_prob_left
        jump_prob_right = jump_prob_right./jump_prob_denominator
        jump_prob_self = jump_prob_self./jump_prob_denominator
        jump_prob_left = jump_prob_left./jump_prob_denominator


        # Compute new boundary value at 'a'
        u[1] = 0 * u0[1] + jump_prob_left[1] * u0[1] +
            jump_prob_left[2] * u0[2] +
            jump_prob_self[1] * u0[1] - ν * u0[1] +
            source_term(x[1],source_params,p0)

        # Compute Interior Points
        u[2:end-1] = jump_prob_right[1:end-2] .* u0[1:end-2] +
            jump_prob_left[3:end] .* u0[3:end] +
            jump_prob_self[2:end-1] .* u0[2:end-1] -
            ν * u0[2:end-1] +
            [source_term(xᵢ, source_params,p0) for xᵢ in x[2:end-1]]

        # The 'mass' at site 'j' at the next time step is the mass at 'j-1'
        # times the probability of right plus the mass at 'j+1' times the
        # probability of jumping left plus the mass at site 'j' that self
        # jumps back to site 'j'.

        # Compute new boundary value at b
        u[end] = jump_prob_right[end-1] * u0[end-1] +
            jump_prob_right[end] * u0[end] +
            jump_prob_self[end] * u0[end] -
            ν * u0[end] + source_term(x[end],source_params,p0)

        u0 = u
        # plot!(1:200, u0, label = "u$i")
    end
    # display(plt)
    return u
end

# a,b,D,V₀, ν, source_params, β, p0, n, T
#@benchmark DTRW_reaction_diffusion(1,100,5.0,0.0,0.5,0.5,0.5,1.0,5.0,1,10,100)
# u = DTRW_reaction_diffusion(1,100,5.0,0.1,0.01,[1.0,0.5],2.0,50.5,201,20)
# # Juno.@enter DTRW_reaction_diffusion(1,100,5.0,0.9,0.5,[1.5,0.5],1.0,50.5,200,10)
# plot(1:201, u)
# # 1+100
# argmin(abs.(u))

# u0 = get_initial_conditions(1,100,0.01,5.0, 0.001, [1.0,0.5], 50.5, 2001)
# plot(1:size(u0,1), u0)
# u0[1001]
# u0[1002]
# u0[1000]
