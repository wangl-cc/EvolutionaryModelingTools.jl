using EvolutionaryModelingTools
using EvolutionaryModelingTools: sample
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

    @testset "Scalar" begin
        include("scalar.jl")
    end

    @testset "Sample" begin
        @test sample(1, 1) == CartesianIndex()
        @test sample(fill(1), 1) == CartesianIndex()
        @test sample([1, 2], 2) == CartesianIndex(2)
        @test sample([1 2], 2) == CartesianIndex(1, 2)
        @test sample(reshape(1:4, 2, 2), 1) == CartesianIndex(1, 1)
        @test sample(reshape(1:4, 2, 2), 2) == CartesianIndex(2, 1)
        @test sample(reshape(1:4, 2, 2), 4) == CartesianIndex(1, 2)
        @test sample(reshape(1:4, 2, 2), 7) == CartesianIndex(2, 2)
    end

    @testset "Tools" begin
        include("tools.jl")
    end
end
