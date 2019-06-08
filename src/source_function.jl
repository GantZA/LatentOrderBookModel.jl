"""
Source Term Function
```julia
source_term(x, source_params, p)
```
where x is a discritzed price point, source_params is a vector of parameters for
the source term and p is the last "observed" mid price point

The function returns a single scalar output

"""

function source_term(
    x::Float64,
    source_params::st_params,
    p::Float64)

    λ = source_params.λ
    μ = source_params.μ
    return λ*tanh(μ*(p-x))
end
