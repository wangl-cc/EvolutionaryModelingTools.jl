const PRESERVED_ARGS = (:ind,)

macro args(ex)
    funcname, funcargs = _parsefunc(Val(ex.head), ex.args...)
    func1 = Expr(
        :(=),
        :($funcname($(PRESERVED_ARGS...), args::NamedTuple,)),
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
    if sym in PRESERVED_ARGS
        return sym
    else
        return Expr(:., :args, QuoteNode(sym))
    end
end
