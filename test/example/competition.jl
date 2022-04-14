module Competition

using Test

include("../../example/competition.jl")

g_result = run_gillespie(MersenneTwister(1), (0.0, T))[1]
g_result_rec = run_gillespie_record(MersenneTwister(1), T)[1]
m_result = run_manually(MersenneTwister(1), T)
m_result_rec = run_manually_record(MersenneTwister(1), T)

@test g_result.v == g_result_rec.v == m_result.v == m_result_rec.v
@test g_result.m == g_result_rec.m == m_result.m == m_result_rec.m

@test all(map(==, getentries(g_result_rec.v), getentries(m_result_rec.v)))
@test all(map(==, values(getentries(g_result_rec.m)), values(getentries(m_result_rec.m))))

end
