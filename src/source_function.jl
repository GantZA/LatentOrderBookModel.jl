# """
# Source Term Function
# ```julia
# source_term = SourceTerm(λ, μ)
# source_term(x, p)
# ```
# where x is a discritzed price point, source_params is a vector of parameters for
# the source term and p is the last "observed" mid price point
#
# The function returns a single scalar output
#
# """

function (st::SourceTerm)(x::Float64,p::Float64)
    return st.λ*tanh(st.μ*(p-x))
end
