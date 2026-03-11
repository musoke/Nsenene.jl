using Nsenene

using Documenter
using Test

doctest(Nsenene)

@testset "Aqua.jl" begin
    include("aqua.jl")
end
