import Nsenene: densities, density
import Nsenene: CylindricalProfile, SphericalProfile

resol = 17
length = 1.0
nfields = 3
m_ = 1:nfields

@testset for profile_type in [CylindricalProfile, SphericalProfile]
    p = profile_type(resol, length, nfields)

    msize = (ones(Int8, ndims(p.psi) - 1)..., nfields)
    m = reshape(m_, msize)

    p.psi .= 1.0

    rhos = densities(p, m)
    rho = density(p, m)

    @test size(rhos) == size(p.psi)
    @test all((rhos .>= 0.0))
    @test selectdim(rhos, ndims(m), 2) == m[2] / m[1] * selectdim(rhos, ndims(m), 1)

    @test size(rho) == size(p.psi)[1:(end - 1)]
    @test all((rho .>= rhos))
    @test all((rho .== sum(m))) # because psi.=1
end
