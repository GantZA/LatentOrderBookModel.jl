using LinearAlgebra

# SPDE
# ϱₜ = -Vₜ * ϱₓ + D * ϱₓₓ - ν * ϱ(x,t) + λ * tanh(μ * (p(t) - x))

function initial_conditions(ds::DTRWSolver, rdpp::ReactionDiffusionPricePath)
    Δx = (ds.b - ds.a)/(rdpp.n_spatial_points-1)
    V₀ = rand(Normal(0, rdpp.σ))
    A = Tridiagonal(
        (V₀/(2*Δx) + rdpp.D/(Δx^2)) * ones(rdpp.n_spatial_points-1),
        ((-2*rdpp.D)/(Δx^2) - rdpp.ν) * ones(rdpp.n_spatial_points),
        (-V₀/(2*Δx) + rdpp.D/(Δx^2)) * ones(rdpp.n_spatial_points-1))

    A[1,2] = 2*rdpp.D/Δx^2
    A[end, end-1] = 2*rdpp.D/Δx^2

    B = .-[rdpp.source_term(xᵢ, ds.mid_price) for xᵢ in ds.x]

    U = A \ B
    return U
end

function (ds::DTRWSolver)(rdpp::ReactionDiffusionPricePath)

    u0 = initial_conditions(ds, rdpp)
    u = u0[:]
    Δx = (ds.b - ds.a)/(rdpp.n_spatial_points-1)
    Δt = (Δx^2) / (2*rdpp.D)

    # Simulate PDE
    for i = 1:rdpp.τ
        ϵ = rand(Normal(0,1))
        Vₜ = rdpp.α*ϵ
        V = -Vₜ*ds.x/(2*rdpp.D)

        jump_prob_right = vcat(exp.(-rdpp.β*V[2:end]),exp(-rdpp.β*V[end]))
        jump_prob_self = exp.(-rdpp.β*V)
        jump_prob_left = vcat(exp(-rdpp.β*V[1]),exp.(-rdpp.β*V[1:end-1]))
        jump_prob_denominator = jump_prob_right + jump_prob_self +
            jump_prob_left

        # Normalizing the probabilities
        jump_prob_right = jump_prob_right./jump_prob_denominator
        jump_prob_self = jump_prob_self./jump_prob_denominator
        jump_prob_left = jump_prob_left./jump_prob_denominator


        # Compute new boundary value at 'a'
        u[1] = jump_prob_left[1] * u0[1] +
            jump_prob_left[2] * u0[2] +
            jump_prob_self[1] * u0[1] - rdpp.ν * u0[1] +
            rdpp.source_term(ds.x[1], ds.mid_price)

        # Compute Interior Points
        u[2:end-1] = jump_prob_right[1:end-2] .* u0[1:end-2] +
            jump_prob_left[3:end] .* u0[3:end] +
            jump_prob_self[2:end-1] .* u0[2:end-1] -
            rdpp.ν * u0[2:end-1] +
            [rdpp.source_term(xᵢ, ds.mid_price) for xᵢ in ds.x[2:end-1]]

        # The 'mass' at site 'j' at the next time step is the mass at 'j-1'
        # times the probability of right plus the mass at 'j+1' times the
        # probability of jumping left plus the mass at site 'j' that self
        # jumps back to site 'j'.

        # Compute new boundary value at b
        u[end] = jump_prob_right[end-1] * u0[end-1] +
            jump_prob_right[end] * u0[end] +
            jump_prob_self[end] * u0[end] -
            rdpp.ν * u0[end] + rdpp.source_term(ds.x[end], ds.mid_price)

        u0 = u[:]
    end
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
