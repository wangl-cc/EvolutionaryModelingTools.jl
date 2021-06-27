"""
    sample([rng=GLOBAL_RNG], A::AbstractArray, [sum_A=sum(A))])

Sampling a index of given array `A` with weight `A`.
"""
function sample end
# sample for single real number
sample(::AbstractRNG, ::Real) = CartesianIndex(1)
sample(::AbstractRNG, ::Real, ::Real) = CartesianIndex(1)
# sample for real scalar
sample(::AbstractRNG, ::AbstractArray{<:Real,0})= CartesianIndex(1)
sample(::AbstractRNG, ::AbstractArray{<:Real,0}, ::Real)= CartesianIndex(1)
# sample for real vector
sample(rng::AbstractRNG, V::AbstractVector{<:Real}, sum_A::Real=sum(V)) =
    CartesianIndex(_sample(rng, V, sum_A))
sample(rng::AbstractRNG, V::Tuple, sum_A::Real=sum(V)) = _sample(rng, V, sum_A)
# sample for real multi-dimensional array
function sample(rng::AbstractRNG, A::AbstractArray{<:Real,N}, sum_A::Real=sum(A)) where {N}
    ind = _sample(rng, vec(A), sum_A)
    return CartesianIndex(ind2sub(A, ind))::CartesianIndex{N}
end

"""
    ind2sub(A::AbstractArray, ind::Integer)

Convert linear index `ind` of  `A` to cartesian index.
"""
ind2sub(A::AbstractArray, ind::Integer) = ind2sub(size(A), ind)
if isdefined(Base, :_ind2sub)
    ind2sub(dims::NTuple, ind::Integer) = Base._ind2sub(dims, ind)
else
    ind2sub(::NTuple{1}, ind::Integer) = ind
    function ind2sub(dims::NTuple, ind::Integer)
        indnext, sub = divrem(ind - 1, dims[1])
        return sub + 1, ind2sub(Base.tail(dims), indnext + 1)...
    end
end

# This method is a modification of `sample([rng], wv::AbstractWeights)` of
# `StatsBase.jl` [MIT License](https://github.com/JuliaStats/StatsBase.jl)
function _sample(rng::AbstractRNG, V, sum_V::Real)
    t = rand(rng) * sum_V
    n = length(V)
    i = 1
    cw = V[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += V[i]
    end
    return i
end
