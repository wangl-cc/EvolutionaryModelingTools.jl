module SIR

using RecordedArrays
using Random
using Test
using EvolutionaryModelingTools

# parameters
clock = ContinuousClock(300.0)
β = 0.001
γ = 0.01

# state
S = recorded(DynamicEntry, clock, 100)
I = recorded(DynamicEntry, clock, 1)
R = recorded(DynamicEntry, clock, 0)

# simulate with this package
@reaction infection begin
    β * S * I
    (S[] -= 1; I[] += 1)
end

@reaction recovery begin
    γ * I
    (I[] -= 1; R[] += 1)
end

ps = (; β, γ, S, I, R)
rs = (infection, recovery)

Random.seed!(1)
init!(clock)

ps_updated = gillespie(clock, ps, rs)[1]

# simulate manually
Random.seed!(1)
init!(clock)

for _ in clock
    infection_rate = β * S * I
    recovery_rate = γ * I

    summed = infection_rate + recovery_rate

    τ = -log(rand()) / summed # compute time interval

    increase!(clock, τ) # update current time

    # sample a reaction and adjust population size
    rn = rand() * summed
    if rn < infection_rate
        S[] -= 1
        I[] += 1
    else
        I[] -= 1
        R[] += 1
    end
end

for sym in (:S, :I, :R)
    @eval @test state($sym) == state(ps_updated.$sym)
    @eval @test getts(getentries($sym)) == getts(getentries(ps_updated.$sym))
    @eval @test getvs(getentries($sym)) == getvs(getentries(ps_updated.$sym))
end

end
