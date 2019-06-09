using LinearAlgebra

include("source_function.jl")

# SPDE
# ϱₜ = -Vₜ * ϱₓ + D * ϱₓₓ - η * ϱ(x,t) + λ * tanh(μ * (p(t) - x))

mutable struct DTRW_RD_Params
    a_b::Array{Float64,1}
    x::Array{Float64,1}
    advective_coeff::Float64
    mid_price::Float64
end


function stochastic_drift(α₀, α, ε)
    return α₀ + α*ε
end

function get_initial_conditions(dtrw_params, rdp_params, st_params)
    a,b         = dtrw_params.a_b
    x           = dtrw_params.x
    n           = rdp_params.n_spatial_points
    V₀          = dtrw_params.advective_coeff
    D           = rdp_params.D
    mid_price   = dtrw_params.mid_price
    η           = rdp_params.η

    Δx = (b - a)/(n-1)
    A = Tridiagonal(
        (V₀/(2*Δx) + D/(Δx^2)) * ones(n-1),
        ((-2*D)/(Δx^2) - η) * ones(n),
        (-V₀/(2*Δx) + D/(Δx^2)) * ones(n-1))

    A[1,2] = 2*D/Δx^2
    A[end, end-1] = 2*D/Δx^2

    B = .-[source_term(xᵢ, st_params, mid_price) for xᵢ in x]

    U = A \ B
    return U
end

function DTRW_reaction_diffusion( dtrw_params, rdp_params, st_params)
    #ϱₜ = -Vₜ * ϱₓ + D * ϱₓₓ - η * ϱ(x,t) + λ * tanh(μ * (p(t) - x))

    # get reaction diffusion path parameters from rdp_params
    a,b         = dtrw_params.a_b
    x           = dtrw_params.x
    V₀          = dtrw_params.advective_coeff
    mid_price   = dtrw_params.mid_price
    β           = rdp_params.boltz_const
    n           = rdp_params.n_spatial_points
    D           = rdp_params.D
    η           = rdp_params.η
    τ           = rdp_params.τ

    u0 = get_initial_conditions(dtrw_params, rdp_params, st_params)
    u = u0[:]
    Δx = (b - a)/(n - 1)
    Δt = (Δx^2) / (2*D)

    # Simulate PDE
    for i = 1:τ
        ϵ = rand(Normal(0,1))
        Vₜ = stochastic_drift(0.0,1.0,ϵ)
        V = -Vₜ*x/(2*D*β)

        jump_prob_right = vcat(exp.(-β*V[2:end]),exp(-β*V[end]))
        jump_prob_self = exp.(-β*V)
        jump_prob_left = vcat(exp(-β*V[1]),exp.(-β*V[1:end-1]))
        jump_prob_denominator = jump_prob_right + jump_prob_self +
            jump_prob_left

        # Normalizing the probabilities
        jump_prob_right = jump_prob_right./jump_prob_denominator
        jump_prob_self = jump_prob_self./jump_prob_denominator
        jump_prob_left = jump_prob_left./jump_prob_denominator


        # Compute new boundary value at 'a'
        u[1] = jump_prob_left[1] * u0[1] +
            jump_prob_left[2] * u0[2] +
            jump_prob_self[1] * u0[1] - η * u0[1] +
            source_term(x[1], st_params, mid_price)

        # Compute Interior Points
        u[2:end-1] = jump_prob_right[1:end-2] .* u0[1:end-2] +
            jump_prob_left[3:end] .* u0[3:end] +
            jump_prob_self[2:end-1] .* u0[2:end-1] -
            η * u0[2:end-1] +
            [source_term(xᵢ, st_params, mid_price) for xᵢ in x[2:end-1]]

        # The 'mass' at site 'j' at the next time step is the mass at 'j-1'
        # times the probability of right plus the mass at 'j+1' times the
        # probability of jumping left plus the mass at site 'j' that self
        # jumps back to site 'j'.

        # Compute new boundary value at b
        u[end] = jump_prob_right[end-1] * u0[end-1] +
            jump_prob_right[end] * u0[end] +
            jump_prob_self[end] * u0[end] -
            η * u0[end] + source_term(x[end], st_params, mid_price)

        u0 = u[:]
    end
    return u
end


# a,b,D,V₀, η, source_params, β, p0, n, T
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
