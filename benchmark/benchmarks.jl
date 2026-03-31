using BenchmarkTools
using Nsenene

const SUITE = BenchmarkGroup()

nfields = 1;
length = 10.0;
resol_r = 1024;
resol_z = 2resol_r;
m = ones(1);

pc = CylindricalProfile(resol_r, resol_z, length, length, nfields);
ps = SphericalProfile(resol_r, length, nfields);

pc.psi[:, :, 1] .= 100 * exp.(-10 * (pc.R .^ 2 .+ pc.z .^ 2));
ps.psi[:, :, 1] .= 100 * exp.(-10 * ps.r .^ 2);

SUITE["total_mass"] = BenchmarkGroup()
SUITE["total_mass"]["cylindrical"] = @benchmarkable total_mass(pc, ones(1, 1, nfields))
SUITE["total_mass"]["spherical"] = @benchmarkable total_mass(ps, ones(1, nfields))

SUITE["gravitational_potential"] = BenchmarkGroup()
SUITE["gravitational_potential"]["cylindrical"] = @benchmarkable gravitational_potential(
    pc, ones(1, 1, nfields)
)
SUITE["gravitational_potential"]["spherical"] = @benchmarkable gravitational_potential(
    ps, ones(1, nfields)
)
