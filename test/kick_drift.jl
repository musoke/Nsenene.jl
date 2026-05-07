import Nsenene: total_mass
import Nsenene: gravitational_potential
import Nsenene: kick!, drift!
import Nsenene: compute_explaps

function radius(p::SphericalProfile)
    return p.r
end

function radius(p::CylindricalProfile)
    return @. (p.R^2 + p.z^2)^0.5
end

resol = 1024
length = 1.0
nfields = 3
m_ = 1:nfields
h = 1e-10

@testset for profile_type in [CylindricalProfile, SphericalProfile]
    profile = profile_type(resol, length, nfields)

    msize = (ones(Int8, ndims(profile.psi) - 1)..., nfields)
    m = reshape(m_, msize)

    r = radius(profile)

    # Define psi dist in first field
    selectdim(profile.psi, ndims(m), 1) .= 100 * exp.(-10r .^ 2)

    density_old = density(profile, m)
    kick!(profile, m, h)
    density_new = density(profile, m)

    # Imaginary time kick should decrease density everywhere
    @test_broken all(density_new .<= density_old)

    # Density should remain 0 in fields that were initially 0
    for field in 2:nfields
        @test iszero(selectdim(profile.psi, ndims(m), field))
    end

    if profile isa CylindricalProfile
        continue
    end

    explaps = compute_explaps(profile, m, h)
    density_old = density(profile, m)
    drift!(profile, m, h, explaps)
    density_new = density(profile, m)

    # Imaginary time drift should decrease density everywhere
    @test_broken all(density_new .<= density_old)
    # Density should remain 0 in fields that were initially 0
    for field in 2:nfields
        @test iszero(selectdim(profile.psi, ndims(m), field))
    end
end
