import Nsenene: ground_state

resol = 2^7
nfields = 2
m_ = 1:nfields
M = [50.0, 50.0]
Lambda = zeros(nfields, nfields)

function reshape_m(m, profile_type)
    if profile_type == SphericalProfile
        msize = (1, size(m)[end])
    elseif profile_type == CylindricalProfile
        msize = (1, 1, size(m)[end])
    end

    return reshape(m_, msize)
end

@testset for profile_type in [CylindricalProfile, SphericalProfile]

    # msize = (ones(Int8, ndims(profile.psi) - 1)..., nfields)
    m = reshape_m(m_, profile_type)

    profile = ground_state(m, Lambda, M, profile_type, resol)

    @test profile isa profile_type
    @test M ≈ total_masses(profile, m)
end
