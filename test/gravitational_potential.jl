import Nsenene: total_mass
import Nsenene: SphericalProfile
import Nsenene: gravitational_potential

function radius(p::SphericalProfile)
    return p.r
end

function radius(p::CylindricalProfile)
    return @. (p.R^2 + p.z^2)^0.5
end

resol = 1024
length = 1.0
nfields = 2
m_ = 1:nfields

@testset for profile_type in [CylindricalProfile, SphericalProfile]
    profile = profile_type(resol, length, nfields)

    msize = (ones(Int8, ndims(profile.psi) - 1)..., nfields)
    m = reshape(m_, msize)

    r = radius(profile)

    selectdim(profile.psi, ndims(m), 1) .= 100 * exp.(-10r .^ 2)

    Phi = gravitational_potential(profile, m)

    M = total_mass(profile, m)
    Phi_approx = -M ./ r

    @test size(Phi) == size(r)

    # TODO: This could break if r[end] is not at large radius
    @test Phi[end] ≈ Phi_approx[end] rtol = 1e-2

    @test all(Phi .<= 0.0)
end
