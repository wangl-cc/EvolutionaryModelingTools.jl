using EvolutionaryModelingTools
using EvolutionaryModelingTools: sample, collectargs, _parsefunc
using Test

@testset "EvolutionaryModelingTools" begin
    @testset "Examples" begin
        # Examples are test with hand-code simulations
        @testset "Logistic growth with noise" begin
            include("example/logistic.jl")
        end
        @testset "Competitive Dynamics with noise and mutation" begin
            include("example/competition.jl")
        end
        @testset "SIR model with noise" begin
            include("example/sir.jl")
        end
    end

    @testset "Sample" begin
        @test sample(1, 1) == CartesianIndex()
        @test sample(fill(1), 1) == CartesianIndex()
        @test sample([1, 2], 2) == CartesianIndex(2)
        @test sample([1  2], 2) == CartesianIndex(1, 2)
        @test sample(reshape(1:4, 2, 2), 1) == CartesianIndex(1, 1)
        @test sample(reshape(1:4, 2, 2), 2) == CartesianIndex(2, 1)
        @test sample(reshape(1:4, 2, 2), 4) == CartesianIndex(1, 2)
        @test sample(reshape(1:4, 2, 2), 7) == CartesianIndex(2, 2)
    end

    @testset "Tools" begin
        @testset "parse function" begin
            fn, args = _parsefunc(:(@inline f(rng, ind, v) = v[ind] = rand(rng)), (:rng, :ind))
            @test fn == :f
            @test args == [:rng, :ind, :(args.v)]
        end

        @testset "collectargs" begin
            @test collectargs(
                quote
                    global exp
                    iv = f(arg) # internal variable
                    mpd = map(exp, iv) # global variable
                    A[1] = 1 # setindex
                    A[begin+1:end] # setindex
                    A[(i=1; begin+i:end-i)]
                    i # i is a local variable defined in above expression
                    A[begin+j:end-j] # j is a arguments needed
                    arg.x
                    g(x) = x + 1 # local function
                    (x, y) -> x + y # anonymous function
                    @inbounds iv[ind] # macro call
                    rand(rng)::Float64 # type annotation
                end
            ) == Set([:A, :arg, :j, :ind, :rng])
            @test collectargs(
                quote
                    iv = f(arg) # internal variable
                    mpd = map(exp, iv) # global variable
                    A[1] = 1 # setindex
                    A[begin+1:end] # setindex
                    arg.x
                    g(x) = x + 1 # local function
                    (x, y) -> x + y # anonymous function
                    @inbounds iv[ind] # macro call
                    rand(rng)::Float64 # type annotation
                end
            ) == Set([:A, :arg, :exp, :ind, :rng]) # if exp is not mark as global, it will be collected
            @test collectargs(
                :((x::Int, y) -> x * y)
            ) == Set([:x, :y])
            @test collectargs(
                :((x, y) -> @ein C[i, j] := x[i] * y[j]),
            ) == Set([:x, :y])
        end
    end
end
