macro args(ex)
    funcname, funcargs = _parsefunc(Val(ex.head), ex.args...)
    func1 = Expr(:(=),
        :($funcname(args::NamedTuple)),
        :($funcname($(funcargs...))),
    )
    return esc(Expr(:block, func1, ex))
end

_parsefunc(::Val{:function}, ex, _...) = _parsefunc(Val(ex.head), ex.args...)
_parsefunc(::Val{:(=)}, ex, _...) = _parsefunc(Val(ex.head), ex.args...)
_parsefunc(::Val{:where}, ex, _...) = _parsefunc(Val(ex.head), ex.args...)
_parsefunc(::Val{:call}, name, args...) = name, _genarg.(args)

_genarg(arg::Symbol) = Expr(:., :args, QuoteNode(arg))
_genarg(arg::Expr) = Expr(:., :args, QuoteNode(arg.args[1]::Symbol))
