using Nsenene
import Nsenene:
    density,
    energy,
    energy_density,
    energy_density_gravity,
    energy_density_kq,
    energy_density_interaction,
    integrate,
    radius

resol = 2^9
rmax = 16.0
nfields = 2
Lambda = zeros(nfields, nfields)

m = 1:nfields
mc = reshape(m, (1, 1, nfields))
ms = reshape(m, (1, nfields))

pc = CylindricalProfile(resol, rmax, nfields)
ps = SphericalProfile(resol, rmax, nfields)

pc.psi[:, :, 1] .= exp.(-10Nsenene.radius(pc) .^ 2)
ps.psi[:, 1] .= exp.(-10Nsenene.radius(ps) .^ 2)

pc.psi[:, :, 2] .= exp.(-10Nsenene.radius(pc) .^ 2)
ps.psi[:, 2] .= exp.(-10Nsenene.radius(ps) .^ 2)

Nsenene.normalize_mass!(pc, mc, [50.0, 40.0])
Nsenene.normalize_mass!(ps, ms, [50.0, 40.0])

@testset "$(typeof(p)) energies have correct shapes and signs" for (p, m) in
                                                                   ((pc, mc), (ps, ms))
    rho_shape = size(density(p, m))

    @test size(energy_density_gravity(p, m)) == rho_shape
    @test integrate(energy_density_gravity(p, m), p) < 0.0

    @test size(energy_density_kq(p, m)) == rho_shape
    @test integrate(energy_density_kq(p, m), p) > 0.0

    @test size(energy_density_interaction(p, m, Lambda)) == rho_shape
    @test integrate(energy_density_interaction(p, m, Lambda), p) == 0.0
    @test integrate(energy_density_interaction(p, m, ones(nfields, nfields)), p) > 0.0
    @test integrate(energy_density_interaction(p, m, -ones(nfields, nfields)), p) < 0.0
    @test integrate(energy_density_interaction(p, m, [0.0 1.0; 1.0 0.0]), p) > 0.0

    @test size(energy_density(p, m, Lambda)) == rho_shape
    @test energy(p, m, Lambda) isa Float64
end

@testset "Consistent energies" begin
    @show E_grav_spherical = integrate(energy_density_gravity(ps, ms), ps)
    @show E_grav_cylindrical = integrate(energy_density_gravity(pc, mc), pc)
    @test E_grav_spherical ≈ E_grav_cylindrical rtol = 0.01
    @show (E_grav_spherical - E_grav_cylindrical) / E_grav_spherical

    @show E_kq_spherical = integrate(energy_density_kq(ps, ms), ps)
    @show E_kq_cylindrical = integrate(energy_density_kq(pc, mc), pc)
    @test E_kq_spherical ≈ E_kq_cylindrical rtol = 0.01

    @show E_int_spherical = integrate(
        energy_density_interaction(ps, ms, ones(nfields, nfields)), ps
    )
    @show E_int_cylindrical = integrate(
        energy_density_interaction(pc, mc, ones(nfields, nfields)), pc
    )
    @test E_int_spherical ≈ E_int_cylindrical rtol = 0.01

    @show E_spherical = energy(ps, ms, ones(nfields, nfields))
    @show E_cylindrical = energy(pc, mc, ones(nfields, nfields))
    @test E_spherical ≈ E_cylindrical rtol = 0.01
end
