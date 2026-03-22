import Nsenene: densities, density
import Nsenene: CylindricalProfile, SphericalProfile

resol = 1000
length = 3.0
nfields = 3
m_ = 1:nfields

@testset "CylindricalProfile" begin
    p = CylindricalProfile(resol, length, nfields)

    msize = (ones(Int8, ndims(p.psi) - 1)..., nfields)
    m = reshape(m_, msize)

    Ms_exact = zeros(nfields)

    for i in 1:nfields
        # psi i is a cylinder of density i
        radius = log(1 + i^2)
        height = log(1 + i)
        p.psi[:, :, i] .= (p.R .< radius) .&& (abs.(p.z) .< height / 2)
        Ms_exact[i] = pi * radius^2 * height * m[i]
    end

    Ms = total_masses(p, m)

    for i in 1:nfields
        @test Ms[i] ≈ Ms_exact[i] rtol = 1e-2
    end

    M = total_mass(p, m)

    @test M ≈ sum(Ms) rtol = 1e-5
end

@testset "SphericalProfile" begin
    p = SphericalProfile(resol, length, nfields)

    msize = (ones(Int8, ndims(p.psi) - 1)..., nfields)
    m = reshape(m_, msize)

    Ms_exact = zeros(nfields)

    for i in 1:nfields
        # psi i is a sphere of density i
        radius = log(1 + i^2)
        p.psi[:, i] .= (p.r .< radius)
        Ms_exact[i] = 4 / 3 * pi * radius^3 * m[i]
    end

    Ms = total_masses(p, m)

    for i in 1:nfields
        @test Ms[i] ≈ Ms_exact[i] rtol = 1e-2
    end

    M = total_mass(p, m)

    @test M ≈ sum(Ms) rtol = 1e-5
end
