module Nsenene

export density
export total_mass, total_masses
export gravitational_potential

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

"""
    total_mass(profile, m)

Calculate the total mass contained in `profile`, with particle masses `m`.
Returns a scalar.

# Examples

```jldoctest
import Nsenene: CylindricalProfile, total_mass

n = 2
m = ones(1, 1, n)

p = CylindricalProfile(1000, 5.0, n)

# psi 1 is a cylinder of density 1.0, radius 1.0, height 2.0
p.psi[:, :, 1] .= (p.R .< 1) .&& (abs.(p.z) .< 2.0/2)

M = total_mass(p, m)
@assert isapprox(M, 2pi, rtol=1e-2)

# psi 2 is a cylinder of density 1.0, radius 1.0, height 4.0
p.psi[:, :, 2] .= (p.R .< 1) .&& (abs.(p.z) .< 4.0/2)

M = total_mass(p, m)
@assert isapprox(M, 2pi + 4pi, rtol=1e-2)

0.0

# output
0.0
```
"""
function total_mass(profile, m)
    return sum(total_masses(profile, m))
end

"""
    total_masses(profile, m)

Calculate the mass contained by each field in `profile`, with particle masses `m`.
Returns an array of masses.

# Examples

```jldoctest
import Nsenene: SphericalProfile, total_masses

n = 3
m = ones(1, n)

p = SphericalProfile(1000, 3.0, n)

# Ball of radius 1
p.psi[:, 1] .= p.r .< 1
# Ball of radius 2
p.psi[:, 2] .= p.r .< 2

M = total_masses(p, m)

@assert isapprox(M[1], 4/3 * π * 1^3, rtol=1e-2)
@assert isapprox(M[2], 4/3 * π * 2^3, rtol=1e-2)
M[3]

# output

0.0
```
"""
function total_masses end

"""
    gravitational_potential(profile, m)

Compute the gravitational potential due to the fields in `profile` with particle masses `m`.
"""
function gravitational_potential end

include("cylindrical.jl")
include("spherical.jl")

import .Cylindrical: CylindricalProfile
import .Spherical: SphericalProfile

end
