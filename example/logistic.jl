using RecordedArrays
using Random
using EvolutionaryModelingTools
using EvolutionaryModelingTools.Scalar

# parameters
const T = 10.0 # max time
const r = 0.5   # growth rate
const d = 0.1   # death rate
const K = 100   # carrying capacity
const n = 10    # initial population size

@cfunc growth_c(r, n) = r * n
@ufunc growth_u!(ind, n) = n[ind] += 1

@cfunc death_c(d, n) = d * n
@ufunc death_u!(ind, n) = n[ind] -= 1

@cfunc comp_c(r, K, n) = r * n * n' / K
@ufunc comp_u!(ind::CartesianIndex{0}, n) = n[ind] -= 1

# the ind must be CartesianIndex{0}, this is just for test
@ufunc check_extinction(ind::CartesianIndex{0}, n) = n[ind] == 0 ? :extinct : :finnish

growth = Reaction(growth_c, growth_u!)
death = Reaction(death_c, death_u!)
competition = Reaction(comp_c, comp_u!)

_maybe_scalar(n::Number) = scalar(n)
_maybe_scalar(n) = n

function run_gillespie(rng=Random.GLOBAL_RNG, T=T, r=r, d=d, K=K, n=n)
    ps = (; r, d, K, n=_maybe_scalar(n))
    rs = (growth, death, competition)
    return gillespie(check_extinction, rng, T, ps, rs)
end

function run_gillespie_record(rng=Random.GLOBAL_RNG, T=T, r=r, d=d, K=K, n=n)
    clock = ContinuousClock(T)
    ps = (; r, d, K, n= recorded(DynamicEntry, clock, n))
    rs = (growth, death, competition)
    return gillespie(check_extinction, rng, clock, ps, rs)
end

function run_manually(rng=Random.GLOBAL_RNG, T=T, r=r, d=d, K=K, n=n)
    t = 0
    while t <= T
        # calculate rate
        growth_rate = r * n              # intrinsic growth
        death_rate = d * n               # intrinsic death
        competition_rate = r * n * n / K # resource competition
        summed = growth_rate + death_rate + competition_rate  # sum a_i
        # break if summed is zero
        summed == 0 && break
        # update current time
        t += -log(rand(rng)) / summed
        # sample a reaction and adjust population size
        rn = rand(rng) * summed
        if rn < growth_rate
            n += 1
        else
            n -= 1
        end
    end
    return n
end

function run_manually_record(rng=Random.GLOBAL_RNG, T=T, r=r, d=d, K=K, n=n)
    clock = ContinuousClock(T)
    n = recorded(DynamicEntry, clock, n)
    for _ in clock
        # calculate rate
        growth_rate = r * n              # intrinsic growth
        death_rate = d * n               # intrinsic death
        competition_rate = r * n * n / K # resource competition
        summed = growth_rate + death_rate + competition_rate  # sum a_i
        # break if summed is zero
        summed == 0 && break
        # update current time
        increase!(clock, -log(rand(rng)) / summed)
        # sample a reaction and adjust population size
        rn = rand(rng) * summed
        if rn < growth_rate
            n[] += 1
        else
            n[] -= 1
        end
    end
    return n
end
