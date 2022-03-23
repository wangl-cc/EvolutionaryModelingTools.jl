using EvolutionaryModelingTools.Scalar

x = scalar(1)
y = scalar(1im)
@test x[] == 1
@test x[1] == 1
@test x[1, 1] == 1
@test x[CartesianIndex()] == 1
@test x[CartesianIndex(1)] == 1
@test x[CartesianIndex(1, 1)] == 1
@test (+x)::Int == 1
@test (-x)::Int == -1
@test (x + 1)::Int == 2
@test (x - 0x1)::Int == 0
@test (x * 1.0)::Float64 == 1.0
@test (x / 1.0f0)::Float32 == 1.0
@test convert(typeof(y), x) isa typeof(y)
@test x + x == 2
@test x - y == 1 - im
@test x * y == im
@test x / y == -im
@test x \ y == im
@test x ^ y == 1
@test real(x)::Int == 1
@test imag(x)::Int == 0
@test real(y) == 0
@test imag(y) == 1
@test float(x)::Float64 == 1.0
@test float(y)::ComplexF64 == im
@test exp(x) == exp(1)
@test ℯ^y == exp(im)
x[] += 0x1
@test x[] == 2
x[1] += 1
@test x[1] == 3
x[1, 1] += 1
@test x[1, 1] == 4
y[CartesianIndex()] *= im
@test y[] == -1
y[CartesianIndex(1)] *= im
@test y[1] == -im
y[CartesianIndex(1, 1)] *= im
@test y[1, 1] == 1
@test exp(x) == exp(4)
@test ℯ^y == exp(1)
@test promote_type(typeof(x), Int) == Int
@test promote_type(typeof(x), Float64) == Float64
@test promote_type(typeof(x), typeof(y)) == Complex{Int}
