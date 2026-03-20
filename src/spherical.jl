module spherical

import LinearAlgebra: Tridiagonal

G = 1

struct SphericalProfile
    r::Vector{Float64}
    psi::Matrix{Complex{Float16}}
end

function SphericalProfile(resol::Integer, length::Real, n::Integer)
    r = range(length / (resol + 1), length, resol)  # FIXME: is this the best way to address r=0?
    psi = zeros(Complex{Float64}, resol, n)

    return SphericalProfile(r, psi)
end

function dr_element(profile::SphericalProfile)
    r = profile.r
    return dr = diff(r)[1]
end

function densities(profile::SphericalProfile, m)
    psi2 = abs2.(profile.psi)

    return out = m .* psi2
end

"""
    density(profile::SphericalProfile, m)

Calculate the total mass density of `profile` with particle masses `m`.
"""
function density(profile::SphericalProfile, m)
    out = sum(densities(profile, m); dims=2)
    return reshape(out, :)
end

"""
    total_masses(r, psi, m)

Calculate the mass contained by each field with radial coordinates `r` and particle masses `m`.

# Examples

```jldoctest
import Nsenene.spherical: total_masses

resol = 1000
nfields = 2

m = ones(1, nfields)
r = range(0, 10, resol)
psi = zeros(resol, nfields)

# Ball of radius 1
psi[r.<1, 1] .= 1

M = total_masses(psi, r, m)

@assert isapprox(M[1], 4/3 * π * 1^3, atol=1e-1)
M[2]

# output

0.0
```
"""
function total_masses(psi, r, m)
    n = size(psi, 1)
    resol = size(psi, 2)

    @assert size(psi, 1) === size(r, 1)
    @assert size(psi, 2) === size(m, 2)

    dr = diff(r)[1]
    M = sum(4π * r .^ 2 .* abs2.(psi); dims=1) .* m * dr

    @assert size(m) === size(M)

    return M[1, :]
end

"""
    total_masses(profile::SphericalProfile, m)

Calculate the mass contained by each field with radial coordinates `r` and particle masses `m`.

# Examples

```jldoctest
import Nsenene.spherical: SphericalProfile, total_masses

n = 3
m = ones(1, n)

p = SphericalProfile(1000, 3.0, n)

# Ball of radius 1
p.psi[:, 1] .= p.r .< 1
p.psi[:, 2] .= p.r .< 2

M = total_masses(p, m)

@assert isapprox(M[1], 4/3 * π * 1^3, rtol=1e-2)
@assert isapprox(M[2], 4/3 * π * 2^3, rtol=1e-2)
M[3]

# output

0.0
```
"""
function total_masses(profile::SphericalProfile, m)
    return total_masses(profile.psi, profile.r, m)
end

function total_mass(psi, r, m)
    return sum(total_masses(psi, r, m))
end

function total_mass(profile::SphericalProfile, m)
    return total_mass(profile.psi, profile.r, m)
end

function d2_dr2(resol)
    out = Tridiagonal(ones(resol - 1), -2 * ones(resol), ones(resol - 1))
    out[resol, resol] = -1.0

    return out
end

function gravitational_potential(profile, m)
    dr = dr_element(profile)
    r = profile.r
    resol = size(r, 1)

    u = similar(r)
    rho = density(profile, m)

    D = d2_dr2(resol)

    u .= D \ reshape(r .* rho, resol)

    u *= 4 * pi * G * dr^2

    return Phi = u ./ r
end

end
