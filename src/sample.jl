"""
    sample([rng=GLOBAL_RNG], A::AbstractArray, [sum_A=sum(A))])

Sampling a index of given array `A` with weight `A`.
"""
function sample(
    ::AbstractRNG,
    x::Real,
    ::Real = x,
)
    return 1
end

function sample(
    rng::AbstractRNG,
    V::Union{AbstractVector,Tuple},
    wv_sum::Real = sum(V),
)
    t = rand(rng) * wv_sum
    n = length(V)
    i = 1
    cw = V[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += V[i]
    end
    return i
end

function sample(
    rng::AbstractRNG,
    A::AbstractArray,
    sum_A::Real = sum(A),
)
    ind = sample(rng, vec(A), sum_A)
    return ind2sub(A, ind)
end

function sample(wa::Union{AbstractArray,Tuple,Real}, wa_sum::Real = sum(wa))
    return sample(Random.GLOBAL_RNG, wa, wa_sum)
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

"""
    montecarlo([rng=GLOBAL_RNG], as)

Monte-Carlo step of Gillespie algorithm.
"""
function montecarlo(rng::AbstractRNG, as::Union{Real,AbstractArray}...)
    as_sum = map(sum, as)
    as_sum_sum = sum(as_sum)
    iszero(as_sum_sum) && return nothing
    τ = -log(rand(rng)) / as_sum_sum
    ind = sample(rng, as_sum, as_sum_sum)
    sub = sample(rng, as[ind], as_sum[ind])
    return τ, ind, sub
end

montecarlo(rs::Union{Real,AbstractArray}...) = montecarlo(Random.GLOBAL_RNG, rs...)
