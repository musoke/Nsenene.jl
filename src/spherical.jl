module Spherical

import LinearAlgebra: Tridiagonal

import ..drift!
import ..densities
import ..density
import ..gravitational_potential
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

function drift!(profile::SphericalProfile, m, h)
    psi = profile.psi

    laplace = laplacian(profile)

    psi .*= exp.(-im * h / 2 ./ m .* (laplace * psi))

    return profile
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

function d2_dr2_vanish_r0(resol)
    out = Tridiagonal(ones(resol - 1), -2 * ones(resol), ones(resol - 1))
    # Assume derivand is 0 at begin-1

    # Asymptote at R=R_max
    out[resol, resol] = -1.0

    return out
end

function d2_dr2_vanish_r0(profile::SphericalProfile)
    resol_r = size(profile.r, 1)

    return d2_dr2_vanish_r0(resol_r)
end

function d1_dr1(resol_r)
    out = 0.5 * Tridiagonal(-ones(resol_r - 1), zeros(resol_r), ones(resol_r - 1))

    # Forward difference at r=begin
    out[begin, begin] = -1
    out[begin, begin + 1] = 1

    # Backward difference at r=end
    out[end, end - 1] = -1
    out[end, end] = 1

    return out
end

function d2_dr2(resol)
    out = Tridiagonal(ones(resol - 1), -2 * ones(resol), ones(resol - 1))
    # Assume derivand has f[begin-1] = f[begin]
    out[begin, begin] = -1.0

    # Asymptote at R=R_max
    out[resol, resol] = -1.0

    return out
end

function d2_dr2(profile::SphericalProfile)
    resol_r = size(profile.r, 1)

    return d2_dr2(resol_r)
end

function laplacian(profile::SphericalProfile)
    r = profile.r
    dr = dr_element(profile)
    resol_r = size(r, 1)

    d1 = d1_dr1(resol_r) / dr
    d2 = d2_dr2(resol_r) / dr^2

    return 2 ./ r .* d1 + d2
end

function gravitational_potential(profile::SphericalProfile, m)
    dr = dr_element(profile)
    r = profile.r
    resol = size(r, 1)

    rho = density(profile, m)
    u = similar(rho)

    D = d2_dr2_vanish_r0(resol)

    u .= D \ reshape(r .* rho, resol)

    u *= 4 * pi * G * dr^2

    return Phi = u ./ r
end

end
