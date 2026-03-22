module Nsenene

export density

export CylindricalProfile
export SphericalProfile

"""
    densities(profile, m)

Calculate the mass density of each field in `profile`, with particle masses `m`.
"""
function densities(profile, m)
    psi2 = abs2.(profile.psi)

    return out = m .* psi2
end

"""
    density(profile, m)

Calculate the total mass density of `profile` with particle masses `m`.
"""
function density end

include("cylindrical.jl")
include("spherical.jl")

import .Cylindrical: CylindricalProfile
import .Spherical: SphericalProfile

end
