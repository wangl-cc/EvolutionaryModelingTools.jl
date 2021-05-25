struct Reaction{C,U}
    c::C # function to calculate a
    u::U # function to update state
end

struct Model{C<:AbstractClock, P<:NamedTuple,R<:Reaction}
    c::C  # clock
    ps::P # params
    rs::R # reactions
end

function oneloop!(rng::AbstractRNG, m::Model)
    ret = montecarlo(rng, (r.c)(m.ps) for r in m.rs)
    return update!(m, state)   
end
oneloop!(m::Model) = oneloop!(Random.GLOBAL_RNG, m)

_update!(::Model, ::Nothing) = true
function _update!(m::Model, state)
    τ, ind, sub = state
    increase!(m.c, τ)
    r = m.rs[ind]
    (r.u)(m.ps)
    return false
end

function gillespie!(rng::AbstractRNG, m::Model)
    term_state = "finnish"
    for _ in m.c
        term = oneloop!(rng, m)
        if term
            term_state = "zero"
            break
        end
    end
    return term_state
end
gillespie!(m::Model) = gillespie!(Random.GLOBAL_RNG, m)
