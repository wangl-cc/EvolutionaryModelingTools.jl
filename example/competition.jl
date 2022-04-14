using RecordedArrays
using RecordedArrays.ResizingTools
using Random
using LoopVectorization
using EvolutionaryModelingTools
using EvolutionaryModelingTools: sample

const T = 10.0
const r = 0.5
const K = 100
const μ = 0.01
const v = [10;]
const m = ones(1, 1)

# simulate with this package
## define reactions

# define growth with @reaction_eq, loop vectorization is disabled but inbounds, fastmath and  offset are enable
@reaction_eq growth r * (1 - μ) v[i] → 2v[i] avx=false inbounds=true fastmath=true offset=true

# define mutation with @reaction and @quickloop with given index
@reaction mutation begin
    @quickloop r * μ * v[i]
    begin
        push!(v, 1)
        n = length(v)
        resize!(m, (n, n))
        # rand competition matrix
        m[:, n] = rand(rng, n)
        m[n, 1:n-1] = rand(rng, 1, n - 1)
    end
end

# define competition with @reaction and @quickloop
@reaction competition begin
    @quickloop (r / K) * v[i] * m[i, j] * v[j] turbo # force enable loop vectorization
    begin
        i = ind[1]
        v[i] -= 1
        if v[i] == 0 # check extinction
            deleteat!(v, i)
            resize!(m, (not(i), not(i)))
        end
    end
end

# run functions for test and benchmark
function run_gillespie(rng=Random.GLOBAL_RNG, T=T, r=r, K=K, μ=μ, v=v, m=m)
    ps = (; r, K, μ, v, m = SimpleRDArray(m))
    return gillespie(rng, T, ps, (growth, mutation, competition))
end

function run_gillespie_record(rng=Random.GLOBAL_RNG, T=T, r=r, K=K, μ=μ, v=v, m=m)
    clock = ContinuousClock(T)
    ps = (;
        r,
        K,
        μ,
        v=recorded(DynamicEntry, clock, v),
        m=recorded(StaticEntry, clock, m),
    )
    return gillespie(rng, clock, ps, (growth, mutation, competition))
end

# simulate manually
function run_manually(rng=Random.GLOBAL_RNG, T=T, r=r, K=K, μ=μ, v=v, m=m)
    t = 0.0
    v′ = deepcopy(v)
    m′ = SimpleRDArray(m)
    while t <= T
        growth_rate = r * (1 - μ) * v′
        mutation_rate = r * μ * v′
        c = r / K # baseline competition coefficient
        competition_rate = @. c * v′ * m′ * v′'
        rates_sum = (sum(growth_rate), sum(mutation_rate), sum(competition_rate))
        rates_sum_acc = accumulate(+, rates_sum)
        rates_sum_sum = rates_sum_acc[end]
        t += -log(rand(rng)) / rates_sum_sum
        rn = rand(rng) * rates_sum_sum
        ind = searchsortedfirst(collect(rates_sum_acc), rn)
        if ind == 1
            i = sample(growth_rate, rn)
            v′[i] += 1
        elseif ind == 2
            i = sample(mutation_rate, rn - rates_sum_acc[1])
            push!(v′, 1)
            n = length(v′)
            resize!(m′, (n, n))
            # rand competition matrix
            m′[:, n] = rand(rng, n)
            m′[n, 1:n-1] = rand(rng, 1, n - 1)
        elseif ind == 3
            i = sample(competition_rate, rn - rates_sum_acc[2])[1]
            v′[i] -= 1
            if v′[i] == 0 # check extinction
                deleteat!(v′, i)
                resize!(m′, (not(i), not(i)))
            end
        end
        length(v′) == 0 && break
    end
    return (; v=v′, m=m′)
end
function run_manually_record(rng=Random.GLOBAL_RNG, T=T, r=r, K=K, μ=μ, v=v, m=m)
    clock = ContinuousClock(T)
    v′ = recorded(DynamicEntry, clock, v)
    m′ = recorded(StaticEntry, clock, m)
    for _ in clock
        growth_rate = r * (1 - μ) * v′
        mutation_rate = r * μ * v′
        c = r / K # baseline competition coefficient
        competition_rate = @. c * v′ * m′ * v′'
        rates_sum = (sum(growth_rate), sum(mutation_rate), sum(competition_rate))
        rates_sum_acc = accumulate(+, rates_sum)
        rates_sum_sum = rates_sum_acc[end]
        τ = -log(rand(rng)) / rates_sum_sum
        increase!(clock, τ)
        rn = rand(rng) * rates_sum_sum
        ind = searchsortedfirst(collect(rates_sum_acc), rn)
        if ind == 1
            i = sample(growth_rate, rn)
            v′[i] += 1
        elseif ind == 2
            i = sample(mutation_rate, rn - rates_sum_acc[1])
            push!(v′, 1)
            n = length(v′)
            resize!(m′, (n, n))
            # rand competition matrix
            m′[:, n] = rand(rng, n)
            m′[n, 1:n-1] = rand(rng, 1, n - 1)
        elseif ind == 3
            i = sample(competition_rate, rn - rates_sum_acc[2])[1]
            v′[i] -= 1
            if v′[i] == 0 # check extinction
                deleteat!(v′, i)
                resize!(m′, (not(i), not(i)))
            end
        end
        length(v′) == 0 && break
    end
    return (; v=v′, m=m′)
end
