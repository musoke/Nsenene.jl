using BenchmarkTools
using Nsenene
using Nsenene: kick!, drift!, compute_explaps_imag

const SUITE = BenchmarkGroup()

nfields = 1;
length = 10.0;
resol_r = 1024;
resol_z = 2resol_r;
m = ones(1);
Lambda = zeros(nfields)
h = 1e-10

pc = CylindricalProfile(resol_r, resol_z, length, length, nfields);
ps = SphericalProfile(resol_r, length, nfields);

exp_laps_c = compute_explaps_imag(pc, ones(1, 1, nfields), h)
exp_laps_s = compute_explaps_imag(ps, ones(1, nfields), h)

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

SUITE["kick!"] = BenchmarkGroup()
SUITE["kick!"]["cylindrical"] = @benchmarkable kick!(pc, ones(1, 1, nfields), Lambda, h)
SUITE["kick!"]["spherical"] = @benchmarkable kick!(ps, ones(1, nfields), Lambda, h)

SUITE["drift!"] = BenchmarkGroup()
SUITE["drift!"]["cylindrical"] = @benchmarkable drift!(
    pc, ones(1, 1, nfields), h, exp_laps_c
)
SUITE["drift!"]["spherical"] = @benchmarkable drift!(ps, ones(1, nfields), h, exp_laps_s)
