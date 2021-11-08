"""
    sample(A::AbstractArray{<:Real,N}, rn::Real)

Select a cartesian index of given weight array `A` and random number `rn`,
where `A` should be an `AbstractArray` with element type of `Real`.

# Example

```jldoctest
julia> A = reshape(1:4, 2, 2)
2×2 reshape(::UnitRange{Int64}, 2, 2) with eltype Int64:
 1  3
 2  4

julia> EvolutionaryModelingTools.sample(A, 1)
CartesianIndex(1, 1)

julia> EvolutionaryModelingTools.sample(A, 2)
CartesianIndex(2, 1)

julia> EvolutionaryModelingTools.sample(A, 4)
CartesianIndex(1, 2)

julia> EvolutionaryModelingTools.sample(A, 7)
CartesianIndex(2, 2)
```
"""
sample
# sample for number
sample(::Real, ::Real) = CartesianIndex(1)
# sample for scalar
sample(::AbstractArray{<:Real,0}, ::Real) = CartesianIndex(1)
# sample for vector
sample(V::AbstractVector{<:Real}, rn::Real) = CartesianIndex(_sample(V, rn))
# sample for multi-dimensional array
function sample(A::AbstractArray{<:Real}, rn::Real)
    ind = _sample(vec(A), rn)
    return CartesianIndex(ind2sub(A, ind))
end

# convert linear index `ind` of  `A` to cartesian index `sub`
ind2sub(A::AbstractArray, ind::Integer) = Base._ind2sub(axes(A), ind)

# Select a index in `1:length(V)` with weight `V` and random number `rn`
# This method is a modification of `sample([rng], wv::AbstractWeights)` of
# `StatsBase.jl` [MIT License](https://github.com/JuliaStats/StatsBase.jl)
# Copyright (c) 2012-2016: Dahua Lin, Simon Byrne, Andreas Noack, Douglas Bates,
# John Myles White, Simon Kornblith, and other contributors.
function _sample(V, rn::Real) # sample with a given random number
    n = length(V)
    i = 1
    cw = V[1]
    while cw < rn && i < n
        i += 1
        @inbounds cw += V[i]
    end
    return i
end
