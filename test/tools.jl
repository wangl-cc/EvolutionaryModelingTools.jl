using Test
using EvolutionaryModelingTools
using EvolutionaryModelingTools: parse_func, parse_eq, parse_args_body, collectaxes
using LoopVectorization

@testset "parse function" begin
    fn, args = parse_func(:(@inline f(rng, ind, v) = v[ind] = rand(rng)), (:rng, :ind))
    @test fn == :f
    @test args == [:rng, :ind, :(args.v)]
end
@testset "@quickloop" begin
    A = rand(4)
    B = rand(4)
    C = rand(4, 4)
    D = rand(4, 4)
    # those functions is sensitive to the order of arguments, which is undetermined
    f1 = @quickloop R1[i, j] := A[i] * B[j]
    f2 = @quickloop R2[j, i] := A[i] * B[j]
    f3 = @quickloop A[i, j] # just make a copy
    f4 = @quickloop R3[i, j, k] := A[i, j] * B[i, k]
    @test (f1(A, B) == A .* B' == transpose(f2(A, B)) ||
           f1(A, B) == B .* A' == transpose(f2(A, B)))
    @test C == f3(C)
    E = f4(C, D)
    equality_f4 = true
    for k in 1:4, j in 1:4, i in 1:4
        equality_f4 &= ((E[i, j, k] == C[i, j] * D[i, k]) ||
            E[i, j, k] == D[i, j] * C[i, k])
    end
    @test equality_f4
end
@testset "collectaxes" begin
    let
        inds, inds_ax = collectaxes(:(C[i, j]), :(A[i] * B[j]))
        @test inds == [:i, :j]
        @test inds_ax[:i] == Dict(:A => [1])
        @test inds_ax[:j] == Dict(:B => [1])
    end
    let
        inds, inds_ax = collectaxes(:(E[]), :(A[i] * C[i, j] * B[j]))
        @test inds == [:i, :j]
        @test inds_ax[:i] == Dict(:A => [1], :C => [1])
        @test inds_ax[:j] == Dict(:B => [1], :C => [2])
    end
    let
        inds, inds_ax = collectaxes(:(E[j, i]), :(A[i] * C[i, j] * B[j]))
        @test inds == [:j, :i]
        @test inds_ax[:i] == Dict(:A => [1], :C => [1])
        @test inds_ax[:j] == Dict(:B => [1], :C => [2])
    end
    let
        inds, inds_ax = collectaxes(:(E[i, j, k]), :(c * A[i] * B[j, k]))
        @test inds == [:i, :j, :k]
        @test inds_ax[:i] == Dict(:A => [1])
        @test inds_ax[:j] == Dict(:B => [1])
        @test inds_ax[:k] == Dict(:B => [2])
    end
    let
        inds, inds_ax = collectaxes(:(E[i, j]), :(A[i] * B[j, j]))
        @test inds == [:i, :j]
        @test inds_ax[:i] == Dict(:A => [1])
        @test inds_ax[:j] == Dict(:B => [1, 2])
    end
    @test_throws ArgumentError collectaxes(:(C[i, k]),  :(A[i] * B[j]))
end
@testset "parse_eq" begin
    @test parse_eq(:(S + I[i] --> 2*I[i])) == (
        Dict{Any, Int}(:S => 1, :(I[i]) => 1),
        Dict{Any, Int}(:S => -1, :(I[i]) => 1)
    )
    @test parse_eq(:(S --> ∅)) == (
        Dict{Any, Int}(:S => 1),
        Dict{Any, Int}(:S => -1)
    )
    @test parse_eq(:(0 --> S)) == (
        Dict{Any, Int}(),
        Dict{Any, Int}(:S => 1)
    )
    @test parse_eq(:(S[i] --> S[i] + S[i])) == (
        Dict{Any, Int}(:(S[i]) => 1),
        Dict{Any, Int}(:(S[i]) => 1)
    )
    @test parse_eq(:(S + I --> 2(S + I))) == (
        Dict{Any, Int}(:S => 1, :I => 1),
        Dict{Any, Int}(:S => 1, :I => 1)
    )
    @test_throws ArgumentError parse_eq(:(S -> S))
    @test_throws ArgumentError("unexpected expression: f(S)") parse_eq(:(S --> f(S)))
end
@testset "parse_args_body" begin
    @test parse_args_body(quote
        global exp
        iv = f1(arg) # internal variable
        mpd = map(exp, iv) # global variable
        A[1] = 1 # setindex
        A[begin+1:end] # getindex
        A[(i = 1; begin+i:end-i)]
        i # i is a local variable defined in above expression
        A[begin+j:end-j] # j is a arguments needed
        arg.x
        f2(x) = x + 1 # local function
        (x, y) -> x + y # anonymous function
        @inbounds iv[ind] # macro call
        rand(rng)::Float64 # type annotation
    end)[1] == Set([:A, :arg, :j, :ind, :rng])
    @test parse_args_body(quote
        iv = f1(arg) # internal variable
        mpd = map(exp, iv) # global variable
        A[1] = 1 # setindex
        A[begin+1:end] # getindex
        arg.x
        f2(x) = x + 1 # local function
        (x, y) -> x + y # anonymous function
        @inbounds iv[ind] # macro call
        rand(rng)::Float64 # type annotation
    end)[1] == Set([:A, :arg, :exp, :ind, :rng]) # if exp is not mark as global, it will be collected
    @test parse_args_body(quote
        x = x.v # redefine argument
        f2(x) = x + y # local argument x and argument y
        (x, y) -> x + y + z # anonymous function use argument z
    end)[1] == Set([:x, :y, :z])
    @test parse_args_body(:((x::Int, y) -> x * y)) == (Set([:x, :y]), :(x * y))
    @test parse_args_body(:((x, y) -> @ein C[i, j] := x[i] * y[j]))[1] == Set([:x, :y])
end
