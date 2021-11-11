module Competition

using RecordedArrays
using Random
using Test
using EvolutionaryModelingTools
using EvolutionaryModelingTools: sample

# parameters
clock = ContinuousClock(10.0)
r = 0.5   # growth rate
K = 100   # carrying capacity
μ = 0.01  # mutation rate

# state
v = recorded(DynamicEntry, clock, [100])
m = recorded(DynamicEntry, clock, ones(1, 1))

# simulate with this package
@reaction growth begin
    r * (1 - μ) * v
    v[ind] += 1
end

@reaction mutation begin
    @. r * μ * v
    begin
        push!(v, 1)
        n = length(v)
        resize!(m, (n, n))
        # rand competition matrix
        m[:, n] = rand(rng, n)
        m[n, 1:n-1] = rand(rng, 1, n-1)
    end
end

@reaction competition begin
    begin
        c = r / K # baseline competition coefficient
        @. c * v * m * v'
    end
    begin
        i, j = ind.I
        v[i] -= 1
        v[j] += 1
        if v[i] == 0 # check extinction
            deleteat!(v, i)
            resize!(m, (not(i), not(i)))
        end
    end
end

@ufunc check_extinction(ind, v) = println(ind, v)

ps = (; r, K, μ, v, m)
rs = (growth, mutation, competition)

Random.seed!(1)

ps_updated = gillespie(clock, ps, rs)[1]

# simulate manually
Random.seed!(1)
init!(clock)

for _ in clock
    growth_rate = r * (1 - μ) * v
    mutation_rate = r * μ * v
    c = r / K # baseline competition coefficient
    competition_rate = @. c * v * m * v'
    rates_sum = (sum(growth_rate), sum(mutation_rate), sum(competition_rate))
    rates_sum_acc = accumulate(+, rates_sum)
    rates_sum_sum = rates_sum_acc[end]
    τ = -log(rand()) / rates_sum_sum
    increase!(clock, τ)
    rn = rand() * rates_sum_sum
    ind = searchsortedfirst(collect(rates_sum_acc), rn)
    if ind == 1
        i = sample(growth_rate, rn)
        v[i] += 1
    elseif ind == 2
        i = sample(mutation_rate, rn - rates_sum_acc[1])
        push!(v, 1)
        n = length(v)
        resize!(m, (n, n))
        # rand competition matrix
        m[:, n] = rand(n)
        m[n, 1:n-1] = rand(1, n-1)
    elseif ind == 3
        i, j = sample(competition_rate, rn - rates_sum_acc[2]).I
        v[i] -= 1
        v[j] += 1
        if v[i] == 0 # check extinction
            deleteat!(v, i)
            resize!(m, (not(i), not(i)))
        end
    end
    length(v) == 0 && break
end

@test ps_updated.v == v
@test ps_updated.m == m

for (ep, em) in zip(getentries(ps.v), getentries(v))
    @test getts(ep) == getts(em)
    @test getvs(ep) == getvs(em)
end

for (ep, em) in zip(getentries(ps.m), getentries(m))
    @test getts(ep) == getts(em)
    @test getvs(ep) == getvs(em)
end

end
