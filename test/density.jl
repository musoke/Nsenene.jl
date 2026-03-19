import Nsenene.spherical: SphericalProfile, densities, density

resol = 17
nfields = 3
m = 1:nfields

spherical_profile = SphericalProfile(resol, 1.0, nfields)

spherical_profile.psi .= ones(nfields, resol)

rhos = densities(spherical_profile, m)
rho = density(spherical_profile, m)

@test size(rhos) == (nfields, resol)
@test rhos[2, :] == 2 * rhos[1, :]

@test size(rho) == (1, resol)
@test rho[1] == sum(m)
