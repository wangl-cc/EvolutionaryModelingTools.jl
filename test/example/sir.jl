module SIR

using Test

include("../../example/sir.jl")

g_result = run_gillespie(MersenneTwister(1), (0.0, T))[1]
g_result_rec = run_gillespie_record(MersenneTwister(1), T)[1]
m_result = run_manually(MersenneTwister(1), T)
m_result_rec = run_manually_record(MersenneTwister(1), T)

for sym in (:S, :I, :R)
    @eval @test g_result.$sym == g_result_rec.$sym == m_result.$sym == m_result_rec.$sym
    @eval @test getts(getentries(g_result_rec.$sym)) == getts(getentries(m_result_rec.$sym))
    @eval @test getvs(getentries(g_result_rec.$sym)) == getvs(getentries(m_result_rec.$sym))
end

end # module SIRTest
