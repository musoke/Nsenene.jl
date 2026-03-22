module Cylindrical

import LinearAlgebra: Tridiagonal

import ..densities
import ..density

G = 1

"""
    struct CylindricalProfile

# Examples

```jldoctest
import Nsenene.Cylindrical: CylindricalProfile

resol_R = 64
resol_z = 128

length_R = 10.0
length_z = 20.0

nfields = 2

p = CylindricalProfile(resol_R, resol_z, length_R, length_z, nfields)
maximum(p.R)

# output

10.0
```
"""
struct CylindricalProfile
    R::Array{Float64,2}
    z::Array{Float64,2}
    psi::Array{Complex{Float16},3}
end

function CylindricalProfile(
    resol_R::Integer, resol_z::Integer, length_R::Real, length_z::Real, nfields::Integer
)
    R = reshape(range(length_R / (resol_R + 1), length_R, resol_R), :, 1)  # FIXME: is this the best way to address R=0?
    z = reshape(range(-length_z / 2, +length_z / 2, resol_z), 1, :)

    psi = zeros(Complex{Float64}, resol_R, resol_z, nfields)

    return CylindricalProfile(R, z, psi)
end

function CylindricalProfile(resol::Integer, length::Real, nfields::Integer)
    return CylindricalProfile(resol, resol, length, length, nfields)
end

function density(profile::CylindricalProfile, m)
    out = sum(densities(profile, m); dims=3)
    return dropdims(out; dims=3)
end

end
