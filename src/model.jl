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

# type stable map for tuple
@generated function gmap(f, tp::Tuple)
    FTs = tp.parameters
    FN  = length(FTs)
    return Expr(:tuple, (:(f(tp[$i])) for i in 1:FN)...)
end

# generate if-elseif branch for each reaction for type stable code
@generated function applyreaction(rng, rs, ps, as, as_sum, ind)
    RTs = rs.parameters
    RN  = length(RTs)
    ex = Expr(:elseif, :(ind==$RN), quote
            (rs[$RN].u)(sample(rng, as[$RN], as_sum[$RN]), ps); nothing
        end)
    for i in RN-1:-1:2
        ex = Expr(:elseif, :(ind==$i), quote
                (rs[$i].u)(sample(rng, as[$i], as_sum[$i]), ps); nothing
            end, ex)
    end
    return Expr(:if, :(ind==1), quote
            (rs[1].u)(sample(rng, as[1], as_sum[1]), ps); nothing
        end, ex)
end

function gillespie!(rng::AbstractRNG, c::AbstractClock, ps::NamedTuple, rs::Tuple)
    term_state = :finnish
    for _ in c
        as = gmap(r -> (r.c)(ps), rs)
        as_sum = gmap(sum, as)
        as_sum_sum = sum(as_sum)
        if iszero(as_sum_sum)
            term_state = :break
            break
        end
        τ = -log(rand(rng)) / as_sum_sum
        increase!(c, τ)
        ind = sample(rng, as_sum, as_sum_sum)::Int
        applyreaction(rng, rs, ps, as, as_sum, ind)
    end
    return term_state
end

function gillespie(rng::AbstractRNG, c::AbstractClock, ps::NamedTuple, rs::Tuple)
    c = deepcopy(c)
    ps_copy = map(p -> _copy_ps(c, p), ps)
    term_state = gillespie!(rng, c, ps_copy, rs)
    return SSAResults(ps_copy, term_state)
end
gillespie(c::AbstractClock, ps::NamedTuple, rs::Tuple) =
    gillespie(Random.GLOBAL_RNG, c, ps, rs)

_copy_ps(::AbstractClock, x) = copy(x)
_copy_ps(c::AbstractClock, x::RecordedArrays.AbstractRArray) =
    setclock(x, c)
