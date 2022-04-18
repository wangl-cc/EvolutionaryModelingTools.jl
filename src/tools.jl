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
    fname, fargs = parse_func(ex, (:t,))
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
    arguments name `rng` is reserved for random number generator,
    and the argument name `ind` is preserved for index of "reaction".
    Don't use those variable names for other usage
"""
macro ufunc(ex)
    fname, fargs = parse_func(ex, (:t, :ind, :rng))
    adapter = :(Base.@propagate_inbounds $fname(rng, t, ind, args::NamedTuple) =
        $fname($(fargs...)))
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
@ufunc Base.@propagate_inbounds growth_u!(ind, x) = begin
    i = ind[1]
    x[i] += 1
end
```
where arguments of functions were collected from given expression automatically.

!!! note
    If there are global variables, this macro may can not collect arguments correctly.
    Especially, for functions which accept function and types as arguments,
    those functions and types may also be collect as an arguments.
    Thus, these variables must be marked as global variables by `global` before use them,
    even functions and types.
    Besides, type annotation expressions like `x::Vector{Int}`,
    types `Vector` and `Int` will not be collected.
    Avoid to defined your reaction with a type arguments for type annotation.

!!! warning
    The expression follow the same name preserve rule as `@cfunc` and `@ufunc`,
    don't use those variable names for other usage.

This macro cannot parse macro expression.
Thus in some cases, arguments of expression with macro may not be collected correctly.
To avoid this, define reaction with anonymous functions may helpfully:
```julia
@reaction reaction begin
    (A, B) -> @einsum C[i, k] := A[i, j] * B[j, k] # where C, i, j, k should not be collected
    begin
        # do something
    end
end
```
in this case, only arguments of anonymous function (`A` and `B`) will be collected.

Besides, there are another macro [`@quickloop`](@ref), which is similar to `@einsum`,
but it don't follow the einstein convenience which just generate nest loops out of the expression.
Because the `@quickloop` generate a function and arguments are provided,
indexing variables and return variable will not be collected.
"""
macro reaction(name::Symbol, block::Expr)
    block = rmlines(block)
    if block.head == :block && length(block.args) == 2
        cargs, cbody = parse_args_body(block.args[1])
        uargs, ubody = parse_args_body(block.args[2])
        return esc(gen_reaction(name, cargs, cbody, uargs, ubody))
    else
        throw(ArgumentError("can't parse expression"))
    end
end

function gen_reaction(name, cargs, cbody, uargs, ubody)
    cname = Symbol("$(name)_c")
    uname = Symbol("$(name)_u!")
    return Expr(
        :block,
        :(Base.@propagate_inbounds $cname(t, args::NamedTuple) =
            $cname($((_warparg(arg, (:t,)) for arg in cargs)...))),
        :($cname($(cargs...)) = $cbody),
        :(Base.@propagate_inbounds $uname(rng, t, ind, args::NamedTuple) =
            $uname($((_warparg(arg, (:t, :ind, :rng)) for arg in uargs)...))),
        :($uname($(uargs...)) = $ubody),
        Expr(:(=), name, :(Reaction($cname, $uname))),
    )
end

# extract name and arguments of given function expression
function parse_func(ex::Expr, preserve)
    head = ex.head
    if head in (:function, :(=), :where) # unpack function expression
        return parse_func(ex.args[1], preserve)
    elseif head == :macrocall # unpack macro expression
        return parse_func(ex.args[end], preserve)
    elseif head == :call # the args of :call expression is the function name and the arguments
        fname = ex.args[1]
        fargs = map(arg -> parse_arg(arg, preserve), ex.args[2:end])
        return fname, fargs
    else
        error("Unknown expr")
    end
end

# unpack argument names
parse_arg(arg::Symbol, preserve) = _warparg(arg, preserve)
parse_arg(arg::Expr, preserve) = _warparg(arg.args[1]::Symbol, preserve)
function _warparg(sym::Symbol, preserve)
    if sym in preserve
        return sym
    else
        return Expr(:., :args, QuoteNode(sym))
    end
end

# collect arguments from given expression
function parse_args_body(ex)
    if ex.head == :(->)
        return unbox_args!(Set{Symbol}(), ex.args[1]), unblock(ex.args[2])
    elseif ex.head == :macrocall && ex.args[1] == Symbol("@quickloop")
        args = filter(!isline, ex.args)
        return _quickloop(args[2:end]...)
    else
        return collectargs!(Set{Symbol}(), Set{Symbol}(), ex), ex
    end
end

function unbox_args!(args, ex::Expr)
    if ex.head == :(::)
        push!(args, ex.args[1])
    elseif ex.head == :tuple
        for arg in ex.args
            unbox_args!(args, arg)
        end
    else
        error("can't parse expression")
    end
    return args
end
unbox_args!(args, s::Symbol) = (push!(args, s), args)

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
    elseif head == :ref # ignore begin and end in reference
        collectargs!(args, vars, ex.args[1])
        foreach(
            arg -> collectargs!(args, TempUnion(vars, (:begin, :end)), arg),
            ex.args[2:end],
        )
    elseif head in (:call, :macrocall) # function and macro will not be treat as a args
        foreach(arg -> collectargs!(args, vars, arg), ex.args[2:end])
    else # in other cases, collect all arguments
        foreach(arg -> collectargs!(args, vars, arg), ex.args)
    end
    return args
end

# union two collections
# but push! to the first collection
struct TempUnion{S,T}
    s::S
    t::T
end
Base.in(item, s::TempUnion) = (item in s.s) || (item in s.t)
Base.push!(s::TempUnion, item) = push!(s.s, item)

"""
@quickloop ex opts...

Generate a function to calculate given expression with nested loops.

# Example
```julia
@quickloop x[i] * m[i, j] * y[j]
```
will generate a function
```julia
(m, y, x)->begin
    ret = Array{...}(undef, size(m, 1), size(m, 2)) # expression for type inference is omitted
    @fastmath @inbounds for j = axes(m, 2), i = axes(m, 1)
        ret[i, j] = x[i] * m[i, j] * y[j]
    end
    return ret
end
```
where the size of return array is determined automatically.
The size can be specified by
```
@quickloop z[j, i] := x[i] * m[i, j] * y[j]
```
which will generate a function
```julia
(m, y, x)->begin
    z = Array{...}(undef, size(m, 2), size(m, 1)) # expression for type inference is omitted
    @fastmath @inbounds for i = axes(m, 1), j = axes(m, 2)
        z[j, i] = x[i] * m[i, j] * y[j]
    end
    return z
end
```
where the return value is named by `z` and the size of return array is
`(size(m, 2), size(m, 1))` instead of `(size(m, 1), size(m, 2))`.

# Options
Options should be specified like `@quickloop ex opt=true/false` or `@quickloop ex opt`,
if the value of option is omitted, it will be treated as `true`.

- `avx` or `turbo`: use `LoopVectorization` to accelerate the loop, default is `false`,
    unless `LoopVectorization` have be loaded.
- `inbounds`: use `@inbounds` to accelerate the loop, default is `true`, this option is ignored
    if `avx` or `turbo` is enabled.
- `fastmath`: use `@fastmath` to accelerate the loop, default is `true`, this option is ignored
    if `avx` or `turbo` is enabled.
- `offset`: if `true`, input arrays will be treat as offset arrays, default is `false`.

!!! note
    The order of function arguments generated by this macro is undetermined.
    Thus, it's recommender use this macro inside of `@reaction`,
    which detect the order of function arguments automatically.
"""
macro quickloop(args...)
    args, fn_body = _quickloop(args...)
    fn_args = Expr(:tuple, args...)
    return esc(Expr(:(->), fn_args, fn_body))
end

function _quickloop(ex, opts...)
    lhs, rhs = _split_eq(ex)
    inds, inds_ax = collectaxes(lhs, rhs)
    fn_args = collectargs!(Set{Symbol}(), Set{Symbol}(inds), rhs)
    fn_body = gen_loop(inds, inds_ax, fn_args, lhs, rhs, opts...)
    return fn_args, fn_body
end

function _split_eq(ex)
    if ex.head == :(=) || ex.head == :(:=)
        return ex.args[1], ex.args[2]
    else
        return Expr(:ref, gensym(:ret)), ex
    end
end

function collectaxes(lhs, rhs)
    # Dict of index symbol and dict of array and dim of the index
    inds_ax = Dict{Symbol,Dict{Symbol,Vector{Int}}}()
    collectaxes!(rhs, inds_ax)
    if length(lhs.args) == 1
        inds = collect(keys(inds_ax))
        max_pos(ind) = maximum(Iterators.flatten(values(inds_ax[ind])))
        sort!(inds; by=max_pos)
        append!(lhs.args, inds)
    else
        inds = lhs.args[2:end]
        for ind in keys(inds_ax)
            if !in(ind, inds)
                throw(ArgumentError(
                    "index symbol `$ind` don't occur in left hand side `$lhs`"
                ))
            end
        end
    end
    return inds, inds_ax
end

function collectaxes!(ex::Expr, inds_ax)
    if ex.head == :ref
        array_name = ex.args[1]
        for i in 2:length(ex.args)
            ind = ex.args[i]
            if !haskey(inds_ax, ind)
                inds_ax[ind] = Dict{Symbol,Vector{Int}}()
            end
            if !haskey(inds_ax[ind], array_name)
                inds_ax[ind][array_name] = Vector{Int}()
            end
            push!(inds_ax[ex.args[i]][array_name], i - 1)
        end
    else
        for arg in ex.args
            collectaxes!(arg, inds_ax)
        end
    end
    return inds_ax
end
collectaxes!(::Any, axs) = axs

function gen_loop(inds, inds_ax, fn_args, lhs, rhs, opts...)
    if isempty(inds) # scalar expression
        return Expr(:return, rhs)
    end
    mode = copy(DEFAULT_LOOP_MODE)
    for ex in opts
        opt, flag = _parse_opt(ex)
        if opt == :turbo
            mode[1] = flag
        elseif opt == :avx
            mode[1] = flag
        elseif opt == :inbounds
            mode[2] = flag
        elseif opt == :fastmath
            mode[3] = flag
        # mode[4] is used to determine of check axes of given array is matched,
        # and it's removed because of the introduction of `ArrayInterface.indices`
        # this mode[4] maybe reused in the future
        elseif opt == :offset
            mode[5] = flag
        else
            error("unknown loop mode: $opt")
        end
    end
    ret_name = lhs.args[1]
    ret_size = Expr[]
    ret_type = _infer_ret_type(inds, inds_ax, fn_args, rhs)
    loop_head = Expr(:block)
    for (i, ind) in enumerate(inds)
        array_name, dims = first(inds_ax[ind])
        dim = first(dims)
        push!(ret_size, :(size($array_name, $dim)))
        array_tuple = Expr(:tuple, ret_name)
        dim_tuple = Expr(:tuple, i)
        for (a, ds) in inds_ax[ind]
            for d in ds
                push!(array_tuple.args, a)
                push!(dim_tuple.args, d)
            end
        end
        ax = Expr(:call, :indices, array_tuple, dim_tuple)
        pushfirst!(loop_head.args, Expr(:(=), ind, ax))
    end
    ret_init = :($ret_name = Array{$ret_type}(undef, $(ret_size...)))
    rhs_offset = if mode[5]
        MacroTools.postwalk(rhs) do ex
            if ex isa Expr && ex.head == :ref
                ex′ = deepcopy(ex)
                for i in 2:length(ex.args)
                    ex′.args[i] = :(axes($(ex.args[1]), $(i-1))[$(ex.args[i])])
                end
                return ex′
            else
                return ex
            end
        end
    else
        rhs
    end
    loop_body = Expr(:block, Expr(:(=), lhs, rhs_offset))
    loop = Expr(:for, loop_head, loop_body)
    if mode[1]
        loop = :(@turbo $loop)
    else
        if mode[2]
            loop = :(@inbounds $loop)
        end
        if mode[3]
            loop = :(@fastmath $loop)
        end
    end
    body = Expr(:block, ret_init, loop, Expr(:return, ret_name))
    return MacroTools.flatten(body)
end

# turbo, inbounds, fastmath, [removed], offset, ...
const DEFAULT_LOOP_MODE = Bool[0, 1, 1, 0, 0, 0, 0, 0]

_parse_opt(opt::Symbol) = opt, true
function _parse_opt(opt::Expr)
    if opt.head == :(=)
        return opt.args[1]::Symbol, opt.args[2]::Bool
    else
        throw(ArgumentError("options should be `opt` or `opt = value::Bool`"))
    end
end

function _infer_ret_type(inds, inds_ax, fn_args, ex)
    ex′ = Expr(:block, ex)
    for ind in inds
        array_name, dims = first(inds_ax[ind])
        dim = first(dims)
        assignment = :($ind = firstindex($array_name, $dim))
        pushfirst!(ex′.args, assignment)
    end
    args = collect(fn_args)
    promote_op_fn = Expr(:(->), Expr(:tuple, args...), ex′)
    args_type = map(arg -> :(typeof($arg)), args)
    ret_type = :(Base.promote_op($promote_op_fn, $(args_type...)))
    return ret_type # expression for inferring the return type
end

"""
    @reaction_eq name param eq [opts...]

Generate a reaction with the given name, parameters and equation.

# Example
```julia
@reaction growth r X[i] --> 2X[i]
```
is equivalent to
```julia
@reaction growth begin
    r * X
    begin
        i = ind[1]
        X[i] += 1
    end
end
```
"""
macro reaction_eq(name::Symbol, param, eq::Expr, opts...)
    lhs_counts, dif_count = parse_eq(eq)
    inds, cargs, cbody = gen_rate(name, param, lhs_counts, opts...)
    update_ex = gen_update(inds, dif_count)
    uargs, ubody = parse_args_body(update_ex)
    return esc(gen_reaction(name, cargs, cbody, uargs, ubody))
end

function parse_eq(eq::Expr)
    if Meta.isexpr(eq, :(-->), 2)
        return _parse_eq(eq.args[1], eq.args[2])
    elseif Meta.isexpr(eq, :call, 3) && eq.args[1] === :(→)
        return _parse_eq(eq.args[2], eq.args[3])
    else
        throw(ArgumentError("unknown expression: $eq, expected `-->` or `→`"))
    end
end
function _parse_eq(lhs, rhs)
    lhs_counts = countsyms(lhs)
    rhs_counts = countsyms(rhs)
    dif_counts = Dict{Any,Int}()
    syms = union(keys(lhs_counts), keys(rhs_counts))
    for sym in syms
        if !haskey(dif_counts, sym)
            dif_counts[sym] = 0
        end
        rhs_count = get(rhs_counts, sym, 0)
        lhs_count = get(lhs_counts, sym, 0)
        dif_counts[sym] = rhs_count - lhs_count
    end
    return lhs_counts, dif_counts
end

function gen_rate(name::Symbol, param, counts, opts...)
    ex = Expr(:call, :*, param)
    for (sym, n) in counts
        if n != 1
            push!(ex.args, :($sym^$n))
        else
            push!(ex.args, sym)
        end
    end
    lhs = Expr(:ref, Symbol(name, :_rate))
    inds, inds_ax = collectaxes(lhs, ex)
    fn_args = collectargs!(Set{Symbol}(), Set{Symbol}(inds), ex)
    fn_body = gen_loop(inds, inds_ax, fn_args, lhs, ex, opts...)
    return inds, fn_args, fn_body
end

function gen_update(inds, counts)
    update = quote
        $(Expr(:tuple, inds...)) = ind.I
    end
    for (sym, n) in counts
        if n != 0
            push!(update.args, :($(_maybe_ref(sym)) += $n))
        end
    end
    return update
end

_maybe_ref(sym::Symbol) = Expr(:ref, sym)
_maybe_ref(sym::Expr) = sym

function countsyms(ex)
    count_dict = Dict{Any, Int}()
    return countsyms!(count_dict, ex, 1)
end

function countsyms!(count_dict, ex::Expr, n::Int)
    if ex.head === :ref
        if haskey(count_dict, ex)
            count_dict[ex] += n
        else
            count_dict[ex] = n
        end
    elseif ex.head === :call && ex.args[1] == :+
        for arg in ex.args[2:end]
            countsyms!(count_dict, arg, n)
        end
    elseif ex.head === :call && ex.args[1] == :* && length(ex.args) == 3 && ex.args[2] isa Int
        countsyms!(count_dict, ex.args[3], ex.args[2] * n)
    else
        throw(ArgumentError("unexpected expression: $ex"))
    end
    return count_dict
end
function countsyms!(count_dict, sym::Symbol, n::Int)
    if haskey(count_dict, sym)
        count_dict[sym] += n
    elseif sym != :∅
        count_dict[sym] = n
    end
    return count_dict
end
countsyms!(count_dict, ::Any, ::Int) = count_dict
