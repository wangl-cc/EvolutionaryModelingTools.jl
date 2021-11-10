# TODO: parse expressions with MacroTools

"""
    @cfunc ex

Define a "calculation" function with an "adapter" methods used to parse args from model.
"calculate" functions take arguments from system state
and calculate "rate"s determining the probability weight of the reaction be selected.

# Example

For function definition:
```julia
@cfunc @inline growth_c(r, x::Vector) = r * x # function to calculate "growth rate"
```
This macro creates two methods, an "adapter" method
```julia
growth_c(args::NamedTuple) = growth_c(args.r, args.x)
```
and the origin method
```julia
@inline growth_c(r, x::Vector) = r * x
```

!!! info
    For function definition with other macros, put those macros after this macro.

!!! warning
    The argument name `t` is reserved for time.
    Don't use those variable names for other usage
"""
macro cfunc(ex)
    fname, fargs = _parsefunc(ex, (:t,))
    adapter = :(Base.@propagate_inbounds $fname(t, args::NamedTuple) = $fname($(fargs...)))
    return esc(Expr(:block, adapter, ex))
end

"""
    @ufunc ex

Define a "update" function with an "adapter" methods used to parse args from model.
"update" functions take arguments from system state and update system state.

# Example

For function definition:
```julia
@ufunc Base.@propagate_inbounds growth_u!(ind::CartesianIndex{1}, x) = x[ind] += 1
```
This macro creates two methods, an "adapter" method
```julia
growth_u!(ind, args::NamedTuple) = growth_u(ind, args.x)
```
and the origin method
```julia
Base.@propagate_inbounds growth_u!(ind::CartesianIndex{1}, x) = x[ind] += 1
```

!!! info
    For function definition with other macros, put those macros after this macro.

!!! warning
    The argument name `t` is reserved for time,
    and the argument name `ind` is preserved for index of "reaction".
    Don't use those variable names for other usage
"""
macro ufunc(ex)
    fname, fargs = _parsefunc(ex, (:t, :ind))
    adapter = :(
        Base.@propagate_inbounds $fname(t, ind, args::NamedTuple) =
            $fname($(fargs...))
    )
    return esc(Expr(:block, adapter, ex))
end

"""
    @reaction name ex

Define a [`Reaction`](@ref) with the given `name` and `ex`.

```julia
@reaction growth begin
    r * x # "calculate" expression
    begin
        i = ind[1]
        x[i] += 1 # "update" expression
    end
end
```
will create a `Reaction` with a "calculation" function:
```julia
@cfunc Base.@propagate_inbounds growth_c(r, x) = r * x
```
and an "update" function:
```julia
@cfunc Base.@propagate_inbounds growth_u!(ind, x) = begin
    i = ind[1]
    x[i] += 1
end
```
where arguments of functions were collected from given expression automatically.

!!! note
    If there are global variables, this macro may can not collect arguments correctly
    Especially, for functions which accept function and types as arguments,
    those functions and types may also be collect as an arguments.
    Thus, these variables must be marked as global variables by `global` before use them,
    even functions and types.
    Besides, type annotation expressions like `x::Vector{Int}`,
    types `Vector` and `Int` will not be collected.
    Avoid to defined your reaction with a type arguments for type annotation.

!!! warning
    The argument name `ind` is preserved for index of "reaction".
    Don't use it as your variable name for other usage in "update" expression.
"""
macro reaction(name::Symbol, block::Expr)
    if block.head == :block && length(block.args) == 4
        cname = Symbol("$(name)_c")
        uname = Symbol("$(name)_u!")
        cbody = block.args[2]
        ubody = block.args[4]
        cargs = collectargs(cbody)
        uargs = collectargs(ubody)
        return esc(Expr(:block,
            :(Base.@propagate_inbounds $cname(t, args::NamedTuple) =
                $cname($((_warparg(arg, (:t,)) for arg in cargs)...))),
            :(Base.@propagate_inbounds $cname($(cargs...)) = $cbody),
            :(Base.@propagate_inbounds $uname(t, ind, args::NamedTuple) =
                $uname($((_warparg(arg, (:t, :ind)) for arg in uargs)...))),
            :(Base.@propagate_inbounds $uname($(uargs...)) = $ubody),
            Expr(:(=), name, :(Reaction($cname, $uname))),
        ))
    else
        throw(ArgumentError("can't parse expression"))
    end
end

# extract name and arguments of given function expression
function _parsefunc(ex::Expr, preserve)
    head = ex.head
    if head in (:function, :(=), :where) # unpack function expression
        return _parsefunc(ex.args[1], preserve)
    elseif head == :macrocall # unpack macro expression
        return _parsefunc(ex.args[end], preserve)
    elseif head == :call # the args of :call expression is the function name and the arguments
        fname = ex.args[1]
        fargs = map(arg -> _parsearg(arg, preserve), ex.args[2:end])
        return fname, fargs
    else
        error("Unknown expr")
    end
end

# unpack argument names
_parsearg(arg::Symbol, preserve) = _warparg(arg, preserve)
_parsearg(arg::Expr, preserve) = _warparg(arg.args[1]::Symbol, preserve)
function _warparg(sym::Symbol, preserve)
    if sym in preserve
        return sym
    else
        return Expr(:., :args, QuoteNode(sym))
    end
end

# collect arguments from given expression
collectargs(ex) = collectargs!(Set{Symbol}(), Set{Symbol}(), ex)

# `args` is a list of symbols will be treat as arguments name
# `vars` is a list of symbols which are names of internal variable or global variables
collectargs!(args, _, _) = args # ignore when third argument is neither a symbol nor a expr
function collectargs!(args, vars, sym::Symbol) # symbols neither in syms or iv will be collected
    if !in(sym, args) && !in(sym, vars) && sym != :(:)
        push!(args, sym)
    end
    return args
end
function collectargs!(args, vars, _ex::Expr)
    ex = shortdef(_ex) # force convert function definition to short form
    head = ex.head
    if head == :(=) # assignment or function definition
        arg1 = ex.args[1]
        if (arg1 isa Symbol || (arg1 isa Expr && arg1.head == :tuple)) # assignment to one or multi new symbol
            collectargs!(args, vars, ex.args[2]) # must collect exprs right of equal firstly #2
            collectargs!(vars, args, ex.args[1]) # collect internal variables
            return args
        elseif arg1.head == :call # local function definition as short form
            push!(vars, arg1.args[1]) # local function name will be internal variable
            func_args = Set{Symbol}() # local function arguments
            collectargs!(func_args, Set{Symbol}(), arg1) # collect local function arguments
            collectargs!(args, union(vars, func_args), ex.args[2]) # collect variables used in local function
            return args
        else
            collectargs!(args, vars, ex.args[1]) # collect left side of equal
            collectargs!(args, vars, ex.args[2]) # collect right side of equal
            return args
        end
    elseif head == :-> # anonymous function
        func_args = Set{Symbol}() # anonymous function arguments
        collectargs!(func_args, Set{Symbol}(), ex.args[1]) # collect function arguments
        collectargs!(args, union(vars, func_args), ex.args[2]) # collect variables used in anonymous function
        return args
    elseif head == :. # for x.y only collect x, y will be ignore
        return collectargs!(args, vars, ex.args[1])
    elseif head == :global # global variables
        return union!(vars, ex.args)
    elseif head == :(::) # type annotation
        return collectargs!(args, vars, ex.args[1]) # collect left side of :: but ignore right side
    elseif head in (:call, :macrocall) # function and macro will not be treat as a args
        foreach(arg -> collectargs!(args, vars, arg), ex.args[2:end])
    else # in other cases, collect all arguments
        foreach(arg -> collectargs!(args, vars, arg), ex.args)
    end
    return args
end
