# EvolutionaryModelingTools

[![Build Status](https://github.com/wangl-cc/EvolutionaryModelingTools.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/wangl-cc/EvolutionaryModelingTools.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/wangl-cc/EvolutionaryModelingTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/wangl-cc/EvolutionaryModelingTools.jl)
[![GitHub](https://img.shields.io/github/license/wangl-cc/EvolutionaryModelingTools.jl)](https://github.com/wangl-cc/EvolutionaryModelingTools.jl/blob/master/LICENSE)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wangl-cc.github.io/EvolutionaryModelingTools.jl/stable/)
[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wangl-cc.github.io/EvolutionaryModelingTools.jl/dev/)

A simple package provides an easy way to build evolutionary biology models
and simulate them by Gillespie's direct method algorithm.

## Why?

[`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl)
is a brilliant suite solving deferential equations,
including simulating jump process with Gillespie algorithm, and it is a very good choice.
However, it is not suitable for solving differential equations
with "variable length" state,
which is the main reason why I created this package.

For example, in a [SIR model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model),
the host population is composed of "Susceptible", "Infected" and "Recovered".

In a normal case, there are only one type of "host" and one type "virus" in the system.
Thus the state of host population can be represent as a vector `u = [S, I, R]`,
In a complex case, the host population can be composed of many types of "host" and infected by many types of "virus".
In a system with `n` types of "hosts" and `m` types of "viruses",
the state of host population can also be represented as a vector
by concatenating the state of each component of host `u = vcat(S, vec(I), vec(R))`,
where `S` is a vector of length `n` and `I`, `R` are matrixes of size `n × m`.
However, in evolutionary biology,
the "mutation" and "extinction" will change the types of hosts and viruses,
which means the `n` and `m` changes during the evolution,
and the length of the state vector `u` will also change.

## How to use?

The package [`Catalyst.jl`](https://github.com/SciML/Catalyst.jl) provides a simple way
to build biochemical reaction for `DifferentialEquations.jl`.
Similarly, this package provides a macro `@reaction_eq` which generate reaction(s),
with given equation.

For example, the infect reaction of above SIR model with multi-type of viruses can be defined as:
```julia
@reaction_eq infect β S + I[i] --> 2I[i]
```
where `i` donates the virus type.  The equation `S + I[i] --> 2I[i]` means
the an host infected by the virus `i` infect a susceptible host with rate `β`,
then convert the the susceptible host to a infectious host.
This expression will not only generate one reaction but a group of reactions trough the index `i`.

However, the mutation and extinction can not be defined easily by the macro `@reaction_eq` currently.
Thus the alternative macro `@reaction` provides a low-level way to build reaction(s)
```julia
@reaction mutation begin
    @quickloop μ * I[i]
    begin
        i = ind[1]
        I[i] -= 1 # the host individual is converted to another type
        push!(I, 1) # add a new type of infectious host to the system
        push!(α, randn() + α[i]) # the virulence of the new host type is generated randomly with mean `α[i]`
    end
end
```
This expression defined mutation of virus which contains two parts:
1. The `@quickloop μ * I[i]` defines how to calculate the rate of mutation,
2. The `begin ... end` block defines what happens when the host is mutated,
   where `ind` is a preserved variable which is used to store the index of the mutated host.

Once you have defined all reactions, put them together as a tuple:
```julia
reactions = (infect, mutation, ...)
```
and define the initial state and parameters of the system as a named tuple:
```julia
params = (β = 0.1, μ = 0.001, ..., S = scalar(100), I = [1], R = [0], α = [0.5])
```
where the `scalar` will create a type similar to `Number`,
but it can be update in-place like `Ref` like `S[] += 1`.

Once reaction and parameters is defined, you can use `gillespie` to simulate the system:
```julia
max_time = 100 # the maximum time of simulation
params′, t, term = gillespie(max_time, params, reactions)
```
where the `gillespie` function returns an tuple `t, ps′, term`
where `t` is the time when the simulation ends, `params′` is an updated `params`,
and `term` is a flag indicating whether the simulation is finished after the maximum time
or break with given flag.

**Note**: Changes of the state will not be recorded by default,
but you can use my another package `RecordedArrays` to record them, like this:

```julia
using RecordedArrays
c = DiscreteClock(max_time) # clock store information of max_time
S = recorded(DynamicEntry, c, 100) # create a recorded number as S with the clock c
I = recorded(DynamicEntry, c, [1]) # create a recorded vector as I with the clock c
R = recorded(DynamicEntry, c, [0]) # create a recorded vector as R with the clock c
α = recorded(StaticEntry, c, [0.5]) # create a recorded vector as α with the clock c
params = (; β = 0.1, μ = 0.001, ..., S, I, R, α) # create new params with recorded S, I, R, α
gillespie(c, params, reactions) # run the simulation with the clock and new params
```

More information about `RecordedArrays`, see its
[documentation](https://wangl-cc.github.io/RecordedArrays.jl/dev).
