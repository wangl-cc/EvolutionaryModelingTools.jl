module EvolutionaryModelingTools

using RecordedArrays
using RecordedArrays.ArrayInterface: indices
using Random
using MacroTools
using Requires

export @cfunc, @ufunc, @reaction, @quickloop, @reaction_eq
export Reaction, gillespie, gillespie!, indices

include("tools.jl")
include("sample.jl")
include("model.jl")
include("scalar.jl")

function __init__()
    @require LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890" begin
        DEFAULT_LOOP_MODE[1] = true
    end
end

end # module
