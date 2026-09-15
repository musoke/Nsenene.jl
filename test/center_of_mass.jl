import Nsenene: CylindricalProfile
import Nsenene.Cylindrical: centers_of_mass_z

resol = 128
length = 3.0
nfields = 3
m_ = 1:nfields

p = CylindricalProfile(resol, length, nfields)

msize = (ones(Int8, ndims(p.psi) - 1)..., nfields)
m = reshape(m_, msize)

for i in 1:nfields
    # disc of mass at p.z[i * 10]
    p.psi[i * 10, :, i] .= randn(resol)
end

coms = centers_of_mass_z(p, m)
@test size(coms) == (nfields,)

for i in 1:nfields
    @test coms[i] ≈ p.z[i * 10]
end
