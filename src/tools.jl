const PRESERVED_ARGS = (:ind,)

macro cfunc(ex)
    funcname, funcargs = _parsefunc(ex)
    func1 = Expr(
        :(=),
        :($funcname(args::NamedTuple)),
        :($funcname($(funcargs...))),
    )
    return esc(Expr(:block, func1, ex))
end

macro ufunc(ex)
    funcname, funcargs = _parsefunc(ex)
    func1 = Expr(
        :(=),
        :($funcname($(PRESERVED_ARGS...), args::NamedTuple)),
        :($funcname($(funcargs...))),
    )
    return esc(Expr(:block, func1, ex))
end

macro reaction(name::Symbol, block::Expr)
    if block.head == :block && length(block.args) == 4
        cname = Symbol("$(name)_c")
        uname = Symbol("$(name)_u!")
        cbody = block.args[2]
        ubody = block.args[4]
        cargs = collectargs(cbody)
        uargs = collectargs(ubody)
        cfunc = _gencfunc(cname, block.args[2])
        ufunc = _genufunc(uname, block.args[4])
        return esc(Expr(:block,
            Expr(:(=), :($cname(args::NamedTuple)), :($cname($(_genarg.(cargs)...)))),
            Expr(:(=), :($cname($(cargs...))), cbody),
            Expr(:(=), :($uname($(PRESERVED_ARGS...), args::NamedTuple)), :($uname($(_genarg.(uargs)...)))),
            Expr(:(=), :($uname($(uargs...))), ubody),
            Expr(:(=), name, :(Reaction($cname, $uname))),
           ))
    else
        throw(ArgumentError("block"))
    end
end

function _gencfunc(name::Symbol, ex)
    funcargs = collectargs(ex)
    return Expr(
        :(=),
        :($name(args::NamedTuple)),
        ex
    )
end

function _genufunc(name::Symbol, ex)
    funcargs = collectargs(ex)
    Expr(:(=), :($name($)), ex)
end

function _parsefunc(ex::Expr)
    head = ex.head
    if head in (:function, :(=), :where)
        return _parsefunc(ex.args[1])
    elseif head == :call
        funcname, args... = ex.args # for julia v1.6+
        funcargs = _parsearg.(args)
        return funcname, funcargs
    else
        error("Unknown expr")
    end
end

_parsearg(arg::Symbol) = _warparg(arg)
_parsearg(arg::Expr) = _warparg(arg.args[1])

function _warparg(sym::Symbol)
    if sym in PRESERVED_ARGS
        return sym
    else
        return Expr(:., :args, QuoteNode(sym))
    end
end

collectargs(ex) = collectargs!(Set{Symbol}(), ex)

collectargs!(syms, ::Any) = syms
collectargs!(syms, ex::Symbol) = (ex in syms ? syms : push!(syms, ex))
function collectargs!(syms, ex::Expr)
    if ex.head == :call # Callable will not be treat as a args
        args = ex.args[2:end]
    else
        args = ex.args
    end
    for arg in args
        collectargs!(syms, arg)
    end
    return syms
end
