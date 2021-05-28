macro args(ex)
    funcname, funcargs = _parsefunc(Val(ex.head), ex.args...)
    outerargs = filter(x -> isa(x, Symbol), funcargs) # like ind
    func1 = Expr(
        :(=),
        :($funcname(args::NamedTuple, $(outerargs...))),
        :($funcname($(funcargs...))),
    )
    return esc(Expr(:block, func1, ex))
end

_parsefunc(::Val{:function}, ex, _...) = _parsefunc(Val(ex.head), ex.args...)
_parsefunc(::Val{:(=)}, ex, _...) = _parsefunc(Val(ex.head), ex.args...)
_parsefunc(::Val{:where}, ex, _...) = _parsefunc(Val(ex.head), ex.args...)
_parsefunc(::Val{:call}, name, args...) = name, _genarg.(args)

_genarg(arg::Symbol) = _warparg(arg)
_genarg(arg::Expr) = _warparg(arg.args[1])

function _warparg(sym::Symbol)
    if sym == :ind
        return :ind
    else
        return Expr(:., :args, QuoteNode(sym))
    end
end
