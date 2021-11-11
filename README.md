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
the "mutation" and "extinction" will change the types of "hosts" and "viruses",
which means the `n` and `m` changes during the evolution,
and the length of the state vector `u` will also change.

## How to use?

Similarly to `DifferentialEquations.jl`,
you must define reactions firstly (jump in `DifferentialEquations.jl`),
For example, the infect reaction of above SIR model is defined as:

```julia
@cfunc infect_c(β, S, I) = β * S' * I # function to calculate the rate of infection
@ufunc function infect_u!(ind, S, I) # function to update the state when infection happens
    # ind is the index of the reaction selected with weight calculated by `infect_c`
    S[ind[2]] -= 1 # S_j  -> S_j  - 1
    I[ind] += 1    # I_ij -> I_ij + 1
    return nothing
end
infect = Reaction(infect_c, infect_u!)
```

where `@cfunc` is a macro to help you define a function to calculate the rate of reaction,
`@ufunc` is a macro to help you define a function to update the state when reaction happens.
Unlike `DifferentialEquations.jl`, the "calculate" function returns an array of rates,
which allows variable length state and contains a lot sub-reactions.
Additionally, the "update" function accepts a special argument `ind`,
which tells you which sub-reaction is selected and should always be the first argument.
All of these is designed to work with variable length state.

Besides, there is another macro `@reaction` helps you define a reaction more easily.
The below code is equivalent to the above code:

```julia
@reaction infect begin
    β * S' * I
    begin
        S[ind[2]] -= 1
        I[ind] += 1
    end
end
```

Once you have defined all reactions, put them together as a tuple:

```julia
reactions = (infect, ...)
```

and define the initial state and parameters of the system as a named tuple:

```julia
params = (β = 0.1, ..., S = [10], I = fill(0, 1, 1), R = fill(0, 1, 1))
```

Then, you can use `gillespie` to simulate the system:

```julia
max_time = 100 # the maximum time of simulation
gillespie(max_time, reactions, params) # return a tuple of (term_state, params)
```

where the `gillespie` function returns the `ps` after the simulation
and a flag indicate whether the simulation is finished after the maximum time,
or break before the maximum time.

**Note**: Changes of the state will not be recorded by default,
but you can use my another package `RecordedArrays` to record them, like this:

```julia
using RecordedArrays
c = DiscreteClock(max_time) # clock store information of max_time
S = recorded(DynamicEntry, c, [10]) # create a recorded vector as S
I = recorded(DynamicEntry, c, fill(0, 1, 1)) # create a recorded matrix as I
R = recorded(DynamicEntry, c, fill(0, 1, 1)) # create a recorded matrix as R
params = (β = 0.1, ..., S = S, I = I, R = R) # add new S, I, R to params
gillespie(c, params, reactions) # run the simulation with new params, must replace the max_time to clock
```

More information about `RecordedArrays`, see its
[documentation](https://wangl-cc.github.io/RecordedArrays.jl/dev).

Examples and references of this package can be found in
[documentation](https://wangl-cc.github.io/EvolutionaryModelingTools.jl/dev).
