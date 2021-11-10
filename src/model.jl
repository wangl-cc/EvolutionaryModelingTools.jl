"""
    Reaction{C,U}

Contain a "calculate" function and an "update" function.
The "calculate" function calculates the a "rate" of reaction with the system state,
which determines the probability weight of the reaction be selected,
and the "update" function updates the system state,
when this reaction were selected randomly.
"""
struct Reaction{C,U}
    c::C # function to calculate "rate"
    u::U # function to update state
end

# generated type stable map for tuple,
# Base.map for large Tuple is not type stable
# however, this package must process large Tuples
@generated function gmap(f, tp::Tuple)
    N  = length(tp.parameters)
    return Expr(:tuple, (:(f(tp[$i])) for i in 1:N)...)
end

# apply the given reaction to update system state
# generate if-elseif branch for each reaction for type stable code
@generated function applyreaction(t, rs, ps, as, acc, ind, rn)
    N  = length(rs.parameters)  # number of reactions
    ex = Expr(:elseif, :(ind==$N), quote
                (rs[$N].u)(t, sample(as[$N], rn - acc[$(N-1)]), ps); nothing
            end
        ) # last elseif block with the last reaction
    for i in N-1:-1:2
        ex = Expr(:elseif, :(ind==$i), quote
                    (rs[$i].u)(t, sample(as[$i], rn - acc[$(i-1)]), ps); nothing
                end,
                ex
            ) # append elseif block form last to the second
    end
    return Expr(:if, :(ind==1), quote
                (rs[1].u)(t, sample(as[1], rn), ps); nothing
            end,
            ex
        ) # append if block at first
end


"""
    gillespie!(hook!, rng::AbstractRNG, c::ContinuousClock, ps::NamedTuple, rs::Tuple)

Simulate the system using the Gillespie algorithm with the given parameters,
and return the terminate state `:finnish`, `:break` or any other state returns by the `hook!`.
The clock `c` and parameters `ps` will be updated during the simulation.

# Arguments

- `hook!`: a function with similar arguments to ["update" functions](@ref @ufunc),
   and it's recommended to create `hook!` with [`@ufunc`](@ref) macro;
   unlike "update" functions, `hook` will be called after each reaction
   and should return a terminate state used to terminate the simulation if it is not `:finnish`.
- `rng`: a random number generator for generate random numbers;
- `c`: a clock recording time, which must be the reference clock of recorded variables in `ps`;
- `ps`: a NamedTuple contains state, parameters even args used by `hook!` of the system;
- `rs`: a tuple contains reactions, all parameters required by reactions must be in `ps` with same name.
"""
function gillespie!(hook!, rng::AbstractRNG, c::ContinuousClock, ps::NamedTuple, rs::Tuple)
    term_state = :finnish # terminate state
    for t in c
        as = gmap(r -> (r.c)(t, ps), rs)   # calculate "rate" for each reaction
        as_sum = gmap(sum, as)             # sum of "rate" for each reaction
        as_sum_acc = accumulate(+, as_sum) # accumulate of the sum of "rate" for all reactions
        as_sum_sum = as_sum_acc[end]       # sum of "rate" for all reactions
        if iszero(as_sum_sum)              # if all "rate" are zero
            term_state = :break            # mark terminate state as break
            break                          # break the loop
        end
        τ = -log(rand(rng)) / as_sum_sum # calculate τ
        t′ = increase!(c, τ) # increase clock time
        rn = rand(rng) * as_sum_sum # generate random number, use this rand number twice is OK
        ind = _sample(as_sum, rn)::Int # sample reaction index, use _sample instead of sample to return Int
        applyreaction(t′, rs, ps, as, as_sum_acc, ind, rn) # apply reaction of index ind to update system state
        term_state = hook!(t′, ind, ps) # call hook function
        term_state == :finnish || break # break if state is not :finnish
    end
    return term_state
end

"""
    gillespie([hook!, rng::AbstractRNG,] c, ps::NamedTuple, rs::Tuple)

Simulate the system using the Gillespie algorithm with the given parameters,
and return a tuple of updated `ps` and terminate state.

# Arguments

- `hook!`: a function with similar arguments to ["update" functions](@ref @ufunc)
   and recommended to created with [`@ufunc`](@ref) macro;
   unlike "update" functions, `hook` will be called after each reaction
   and should return a terminate state used to terminate the simulation if it is not `:finnish`.
- `rng`: a random number generator for generate random numbers;
- `c`: a `ContinuousClock`, a end time or a tuple of a begin and a end time;
- `ps`: a NamedTuple contains state, parameters even args used by `hook!` of the system;
- `rs`: a tuple contains reactions, all parameters required by reactions must be in `ps` with same name.
"""
function gillespie(hook!, rng::AbstractRNG, c, ps::NamedTuple, rs::Tuple)
    c′, ps′ = _copy_args(c, ps) # copy args to avoid mutation of ps
    term_state = gillespie!(hook!, rng, c′, ps′, rs)
    return ps′, term_state
end
gillespie(c, ps::NamedTuple, rs::Tuple) =
    gillespie((_...) -> :finnish, Random.GLOBAL_RNG, c, ps, rs)
gillespie(rng::AbstractRNG, c,  ps::NamedTuple, rs::Tuple) =
    gillespie((_...) -> :finnish, rng, c, ps, rs)
gillespie(hook!, c, ps::NamedTuple, rs::Tuple) =
    gillespie(hook!, Random.GLOBAL_RNG, c, ps, rs)

_copy_args(c::ContinuousClock, ps) = deepcopy((c, ps))
_copy_args(c::Real, ps) = ContinuousClock(c), deepcopy(ps)
_copy_args((s, e)::Tuple{Real,Real}, ps) = ContinuousClock(e, s), deepcopy(ps)
