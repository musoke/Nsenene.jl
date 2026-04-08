using Nsenene

using Documenter
using Test

doctest(Nsenene)

@testset "density" begin
    include("density.jl")
end

@testset "total mass" begin
    include("total_mass.jl")
end

@testset "derivatives" begin
    include("derivatives.jl")
end

@testset "gravitational potential" begin
    include("gravitational_potential.jl")
end

@testset "kick" begin
    include("kick.jl")
end

@testset "Aqua.jl" begin
    include("aqua.jl")
end
