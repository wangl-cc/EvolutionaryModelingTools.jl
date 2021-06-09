module GillespieTools

using RecordedArrays
using Random

export @cfunc, @ufunc, @reaction
export Reaction, Model, rsample, rupdate!, gillespie, gillespie!
export rsample

include("tools.jl")
include("sample.jl")
include("model.jl")

end # module
