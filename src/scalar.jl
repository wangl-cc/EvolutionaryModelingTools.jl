module Scalar

export scalar

"""
    scalar(x)

Convert `x` to a ScalarType.
"""
scalar

mutable struct NumberScalar{T<:Number} <: Number
    v::T
end
scalar(v::Number) = NumberScalar(v)
@inline Base.getindex(s::NumberScalar) = s.v
@inline Base.setindex!(s::NumberScalar{T}, v::T) where {T<:Number} = s.v = v

mutable struct RealScalar{T<:Real} <: Real
    v::T
end
scalar(v::Real) = RealScalar(v)
@inline Base.getindex(s::RealScalar) = s.v
@inline Base.setindex!(s::RealScalar{T}, v::T) where {T<:Real} = s.v = v

"""
    ScalarType{T}

ScalarType is a type of number used to represent a scalar value in a model.
Its value is accessible via `getindex` and can be updated via `setindex!`.
"""
const ScalarType{T} = Union{NumberScalar{T}, RealScalar{T}}
const ST = ScalarType

Base.convert(::Type{T}, s::ST{T}) where {T<:Number} = s[]
Base.convert(::Type{T}, s::ST) where {T<:Number} = convert(T, s[])
Base.convert(::Type{T}, s::ST) where {P, T<:ST{P}} = T(convert(P, s[]))

Base.promote_rule(::Type{<:ST{S}}, ::Type{T}) where {S<:Number, T<:Number} =
    promote_type(S, T)
Base.promote_rule(::Type{<:ST{S}}, ::Type{<:ST{T}}) where {S<:Number, T<:Number} =
    promote_type(S, T)

# getindex
@inline function Base.getindex(x::ST, i::Integer)
    @boundscheck isone(i) || throw(BoundsError(x, i))
    return x[]
end
@inline function Base.getindex(x::ST, I::Integer...)
    @boundscheck all(isone, I) || throw(BoundsError(x, I))
    return x[]
end
@inline Base.getindex(x::ST, ::CartesianIndex{0}) = x[] # to fix ambiguity
Base.@propagate_inbounds Base.getindex(x::ST, I::CartesianIndex) = x[I.I...]

# setindex
@inline Base.setindex!(x::ST{T}, v) where {T} = setindex!(x, convert(T, v)) # type conversion
@inline function Base.setindex!(x::ST, v, i::Integer)
    @boundscheck isone(i) || throw(BoundsError(x, i))
    return setindex!(x, v)
end
@inline function Base.setindex!(x::ST, v, I::Integer...)
    @boundscheck all(isone, I) || throw(BoundsError(x, I))
    return setindex!(x, v)
end
@inline Base.setindex!(x::ST, v, ::CartesianIndex{0}) = setindex!(x, v) # to fix ambiguity
Base.@propagate_inbounds Base.setindex!(x::ST, v, I::CartesianIndex) = setindex!(x, v, I.I...)

# mathematical operations
for op in (:+, :-, :conj, :real, :imag, :float, :exp, :log)
    @eval @inline Base.$op(x::ST) = $op(x[])
end

for op in (:+, :-, :*, :/, :\, :^, :(==))
    @eval @inline Base.$op(x::ST{T}, y::ST{T}) where {T} = $op(x[], y[])
end

end # module Scalar
