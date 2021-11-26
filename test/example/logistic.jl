module Logistic

using RecordedArrays
using Random
using Test
using EvolutionaryModelingTools

const c = ContinuousClock(100.0) # define a clock, the population will growth for 100 time unit

const r = 0.5
const d = 0.1
const K = 100
const n = recorded(DynamicEntry, c, 10)   # define a scalar to record population size

# simulate with this package

Random.seed!(1)

@cfunc growth_c(r, n) = r * n
@ufunc growth_u!(ind, n) = n[ind.I...] += 1 # ensure the indexing is correct

@cfunc death_c(d, n) = d * n
@ufunc death_u!(ind, n) = n[ind.I...] -= 1 # ensure the indexing is correct

@cfunc comp_c(r, K, n) = r * n * n / K
@ufunc comp_u!(ind::CartesianIndex{0}, n) = n[ind] -= 1

# the ind must be CartesianIndex{0}, this is just for test
@ufunc check_extinction(ind::CartesianIndex{0}, n) = n[ind] == 0 ? :extinct : :finnish # check

growth = Reaction(growth_c, growth_u!)
death = Reaction(death_c, death_u!)
competition = Reaction(comp_c, comp_u!)

ps = (; r, d, K, n)
rs = (growth, death, competition)

e1 = getentries(gillespie(check_extinction, c, ps, rs)[1].n)

# simulate manually

Random.seed!(1)

for _ in c
    global n
    # evaluate a_i
    growth_rate = r * n              # intrinsic growth
    death_rate = d * n               # intrinsic death
    competition_rate = r * n * n / K # resource competition

    summed = growth_rate + death_rate + competition_rate  # sum a_i

    τ = -log(rand()) / summed # compute time interval

    increase!(c, τ) # update current time

    # sample a reaction and adjust population size
    rn = rand() * summed
    if rn < growth_rate
        n[1] += 1
    else
        n[1] -= 1
    end

    state(n) <= 0 && break # break if population extinct
end

e2 = getentries(n)

@test getts(e1) == getts(e2)
@test getvs(e1) == getvs(e2)

end
