module Spherical

import LinearAlgebra: Tridiagonal

import ..densities
import ..density
import ..total_masses

G = 1

struct SphericalProfile
    r::Vector{Float64}
    psi::Matrix{Complex{Float64}}
end

function SphericalProfile(resol::Integer, length::Real, nfields::Integer)
    r = range(length / (resol + 1), length, resol)  # FIXME: is this the best way to address r=0?
    psi = zeros(Complex{Float64}, resol, nfields)

    return SphericalProfile(r, psi)
end

function dr_element(profile::SphericalProfile)
    r = profile.r
    return dr = diff(r)[1]
end

function density(profile::SphericalProfile, m)
    out = sum(densities(profile, m); dims=2)
    return dropdims(out; dims=2)
end

"""
    _total_masses(r, psi, m)

Calculate the mass contained by each field of `psi` with radial coordinates `r` and particle masses `m`.

# Examples

```jldoctest
import Nsenene.Spherical: _total_masses

resol = 1000
nfields = 2

m = ones(1, nfields)
r = range(0, 10, resol)
psi = zeros(resol, nfields)

# Ball of radius 1
psi[r.<1, 1] .= 1

M = _total_masses(psi, r, m)

@assert isapprox(M[1], 4/3 * π * 1^3, atol=1e-1)
M[2]

# output

0.0
```
"""
function _total_masses(psi, r, m)
    @assert size(psi, 1) === size(r, 1)
    @assert size(psi, 2) === size(m, 2)

    dr = diff(r)[1]
    M = sum(4π * r .^ 2 .* abs2.(psi); dims=1) .* m * dr

    @assert size(m) === size(M)

    return reshape(M, :)
end

function total_masses(profile::SphericalProfile, m)
    return _total_masses(profile.psi, profile.r, m)
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

    rho = density(profile, m)
    u = similar(rho)

    D = d2_dr2(resol)

    u .= D \ reshape(r .* rho, resol)

    u *= 4 * pi * G * dr^2

    return Phi = u ./ r
end

end
