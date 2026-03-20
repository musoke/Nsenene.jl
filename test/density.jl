import Nsenene.spherical: SphericalProfile, densities, density

resol = 17
nfields = 3
m = reshape(1:nfields, 1, nfields)

spherical_profile = SphericalProfile(resol, 1.0, nfields)

spherical_profile.psi .= 1.0

rhos = densities(spherical_profile, m)
rho = density(spherical_profile, m)

@test size(rhos) == (resol, nfields)
@test rhos[:, 2] == 2 * rhos[:, 1]

@test size(rho) == (resol,)
@test rho[1] == sum(m)
