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
    N = length(tp.parameters)
    return Expr(:tuple, (:(f(tp[$i])) for i in 1:N)...)
end

# generated type stable accumulate for tuple
# Base.accumulate for tuple return a tuple
# this method return a Vector
@generated function gaccumulate(tp::Tuple)
    ret_type = promote_type(tp.parameters...)
    N = length(tp.parameters)
    return quote
        ret = Vector{$ret_type}(undef, $N)
        ret[1] = tp[1]
        $((:(@inbounds ret[$i] = ret[$(i - 1)] + tp[$i]) for i in 2:N)...)
        return ret
    end
end

# apply the given reaction to update system state
# generate if-elseif branch for each reaction for type stable code
# as[i], rs[i] are type unstable without @generated
@generated function applyreaction(hook!, rng::AbstractRNG, t, rs, ps, as, acc, ind, rn)
    N = length(rs.parameters)  # number of reactions
    sub = Symbol("sub_$N")
    ex = Expr(
        :elseif,
        :(ind == $N),
        quote
            $sub = sample(as[$N], rn - acc[$(N - 1)])
            (rs[$N].u)(rng, t, $sub, ps)
            return hook!(rng, t, $sub, ps)
        end,
    ) # last elseif block with the last reaction
    for i in N-1:-1:2
        sub = Symbol("sub_$i")
        ex = Expr(
            :elseif,
            :(ind == $i),
            quote
                $sub = sample(as[$i], rn - acc[$(i - 1)])
                (rs[$i].u)(rng, t, $sub, ps)
                return hook!(rng, t, $sub, ps)
            end,
            ex,
        ) # append elseif block form last to the second
    end
    return Expr(:if, :(ind == 1), quote
        sub_1 = sample(as[1], rn)
        (rs[1].u)(rng, t, sub_1, ps)
        return hook!(rng, t, sub_1, ps)
    end, ex) # append if block at first
end

"""
    gillespie!(hook!, rng::AbstractRNG, c::ContinuousClock, ps::NamedTuple, rs::Tuple)

Simulate the system using the Gillespie algorithm with the given parameters,
and return a tuple of end time
and the terminate state `:finnish`, `:zero` or any other state returns by the `hook!`.
The terminate state `:finnish` means that simulation reach to the end time,
and `:zero` means the simulation break because the total "reaction rate" is zero,
besides, `hook!` should return a symbol terminate state like `:break`.
if the return value of `hook!` is not `:finnish`, the simulation will be terminated.
The clock `c` and parameters `ps` will be updated during the simulation.

# Arguments

- `hook!`: a function with similar arguments to "update" functions,
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
    t = currenttime(c)
    for _ in c
        # use t before update it may cause type un-stability
        as = gmap(r -> (r.c)(currenttime(c), ps), rs)   # calculate "rate" for each reaction
        as_sum = gmap(sum, as)             # sum of "rate" for each reaction
        as_sum_acc = gaccumulate(as_sum)   # accumulate of the sum of "rate" for all reactions
        as_sum_sum = as_sum_acc[end]       # sum of "rate" for all reactions
        if iszero(as_sum_sum)              # if all "rate" are zero
            term_state = :zero             # mark terminate state as break
            break                          # break the loop
        end
        τ = -log(rand(rng, typeof(as_sum_sum))) / as_sum_sum   # calculate τ
        t = increase!(c, τ) # increase clock time
        rn = rand(rng, typeof(as_sum_sum)) * as_sum_sum # generate random number, use this rand number twice is OK
        ind = searchsortedfirst(as_sum_acc, rn)::Int # sample reaction index, use _sample instead of sample to return Int
        term_state = applyreaction(hook!, rng, t, rs, ps, as, as_sum_acc, ind, rn) # apply reaction of index ind to update system state
        term_state == :finnish || break # break if state is not :finnish
    end
    return t, term_state
end

"""
    gillespie([hook!, rng::AbstractRNG,] c, ps::NamedTuple, rs::Tuple)

Simulate the system using the Gillespie algorithm with the given parameters,
and return a tuple of updated `ps`, end time and terminate state.
More about terminate state, see [`gillespie!`](@ref).

# Arguments

- `hook!`: a function with similar arguments to "update" functions
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
    t, term_state = gillespie!(hook!, rng, c′, ps′, rs)
    return ps′, t, term_state
end
gillespie(c, ps::NamedTuple, rs::Tuple) =
    gillespie((_...) -> :finnish, Random.GLOBAL_RNG, c, ps, rs)
gillespie(rng::AbstractRNG, c, ps::NamedTuple, rs::Tuple) =
    gillespie((_...) -> :finnish, rng, c, ps, rs)
gillespie(hook!, c, ps::NamedTuple, rs::Tuple) =
    gillespie(hook!, Random.GLOBAL_RNG, c, ps, rs)

_copy_args(c::ContinuousClock, ps) = deepcopy((c, ps))
_copy_args(c::Real, ps) = ContinuousClock(c), deepcopy(ps)
_copy_args((s, e)::Tuple{Real,Real}, ps) = ContinuousClock(e, s), deepcopy(ps)
