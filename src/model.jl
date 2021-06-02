using RecordedArrays: AbstractClock

struct Reaction{C,U}
    c::C # function to calculate a
    u::U # function to update state
end

struct Model{G<:AbstractRNG,C<:AbstractClock,P<:NamedTuple,R<:Tuple}
    rng::G # rng
    c::C   # clock
    ps::P  # params
    rs::R  # reactions
end

struct SSAResults{T}
    ps::T
    state::Symbol
end

function rsample(m::Model)
    rng = m.rng
    as = [(r.c)(m.ps) for r in m.rs]
    as_sum = map(sum, as)
    as_sum_sum = sum(as_sum)
    iszero(as_sum_sum) && return nothing
    τ = -log(rand(rng)) / as_sum_sum
    ind = sample(rng, as_sum, as_sum_sum)
    sub = sample(rng, as[ind], as_sum[ind])
    return τ, ind, sub
end

rupdate!(::Model, ::Nothing) = true
function rupdate!(m::Model, state)
    τ, ind, sub = state
    increase!(m.c, τ)
    r = m.rs[ind]
    (r.u)(sub, m.ps)
    return false
end

function gillespie!(m::Model)
    term_state = :finnish
    for _ in m.c
        term = rupdate!(m, rsample(m))
        if term
            term_state = :break
            break
        end
    end
    return term_state
end

function gillespie(rng::AbstractRNG, c::AbstractClock, ps::NamedTuple, rs::Tuple)
    c = deepcopy(c)
    ps_copy = typeof(ps)(_copy_ps(c, p) for p in ps)
    m = Model(rng, c, ps_copy, rs)
    term_state = gillespie!(m)
    return SSAResults(m.ps, term_state)
end
gillespie(c::AbstractClock, ps::NamedTuple, rs::Tuple) =
    gillespie(Random.GLOBAL_RNG, c, ps, rs)

_copy_ps(::AbstractClock, x) = x
_copy_ps(c::AbstractClock, x::RecordedArrays.AbstractRArray) =
    setclock(x, c)
