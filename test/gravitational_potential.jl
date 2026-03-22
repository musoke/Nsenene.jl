import Nsenene: total_mass
import Nsenene: SphericalProfile
import Nsenene.Spherical: gravitational_potential

resol = 1_000
nfields = 3
m = reshape(1:nfields, 1, nfields)

profile = SphericalProfile(resol, 1.0, nfields)
r = profile.r

profile.psi[:, 1] .= 100 * exp.(-10r .^ 2)

Phi = gravitational_potential(profile, m)

M = total_mass(profile, m)
Phi_approx = -M ./ r

@test size(Phi) == size(profile.r)

@test Phi[end] ≈ Phi_approx[end] rtol = 1e-2
