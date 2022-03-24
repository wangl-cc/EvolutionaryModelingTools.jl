using RecordedArrays
using Random
using EvolutionaryModelingTools
using EvolutionaryModelingTools.Scalar

# parameters
const T = 100.0
const β = 0.001
const γ = 0.01
const S = 100
const I = 1
const R = 0

## simulation with this package
@reaction infection begin
    β * S * I
    (S[] -= 1; I[] += 1)
end

@reaction recovery begin
    γ * I
    (I[] -= 1; R[] += 1)
end

run_gillespie(rng=Random.GLOBAL_RNG, T=T, β=β, γ=γ, S=S, I=I, R=R) = gillespie(
    rng,
    T,
    (; β, γ, S=scalar(S), I=scalar(I), R=scalar(R)),
    (infection, recovery),
)

function run_gillespie_record(rng=Random.GLOBAL_RNG, T=T, β=β, γ=γ, S=S, I=I, R=R)
    clock = ContinuousClock(T)
    ps = (;
        β,
        γ,
        S=recorded(DynamicEntry, clock, S),
        I=recorded(DynamicEntry, clock, I),
        R=recorded(DynamicEntry, clock, R),
    )
    return gillespie(rng, clock, ps, (infection, recovery))
end

## simulation manually
function run_manually(rng=Random.GLOBAL_RNG, T=T, β=β, γ=γ, S=S, I=I, R=R)
    t = 0
    while t <= T
        # calculate rates
        infection_rate = β * S * I
        recovery_rate = γ * I
        summed = infection_rate + recovery_rate
        # break if summed is zero
        summed == 0 && break
        # update current time
        t += -log(rand(rng)) / summed
        # sample a reaction and adjust population size
        rn = rand(rng) * summed
        if rn < infection_rate
            S -= 1
            I += 1
        else
            I -= 1
            R += 1
        end
    end
    return (; t, S, I, R)
end

function run_manually_record(rng=Random.GLOBAL_RNG, T=T, β=β, γ=γ, S=S, I=I, R=R)
    clock = ContinuousClock(T)
    S′ = recorded(DynamicEntry, clock, S)
    I′ = recorded(DynamicEntry, clock, I)
    R′ = recorded(DynamicEntry, clock, R)
    for _ in clock
        # calculate rates
        infection_rate = β * S′ * I′
        recovery_rate = γ * I′
        summed = infection_rate + recovery_rate
        # break if summed is zero
        summed == 0 && break
        # update current time
        increase!(clock, -log(rand(rng)) / summed)
        # sample a reaction and adjust population size
        rn = rand(rng) * summed
        if rn < infection_rate
            S′[] -= 1
            I′[] += 1
        else
            I′[] -= 1
            R′[] += 1
        end
    end
    return (; clock, S=S′, I=I′, R=R′)
end
