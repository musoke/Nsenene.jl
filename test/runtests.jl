using Nsenene

using Documenter
using Test

doctest(Nsenene)

@testset "density" begin
    include("density.jl")
end

@testset "gravitational potential" begin
    include("gravitational_potential.jl")
end

@testset "Aqua.jl" begin
    include("aqua.jl")
end
