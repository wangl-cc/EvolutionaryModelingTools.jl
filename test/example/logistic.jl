module Logistic

using Test

include("../../example/logistic.jl")

g_result = run_gillespie(MersenneTwister(1))[1].n
g_result_record = run_gillespie_record(MersenneTwister(1))[1].n
m_result = run_manually(MersenneTwister(1))
m_result_record = run_manually_record(MersenneTwister(1))

@test g_result == g_result_record == m_result == m_result_record
e1 = getentries(g_result_record)
e2 = getentries(m_result_record)
@test getts(e1) == getts(e2)
@test getvs(e1) == getvs(e2)

end
