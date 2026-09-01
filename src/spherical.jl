module Spherical

import LinearAlgebra: Tridiagonal

import ..compute_explaps_imag
import ..drift!
import ..densities
import ..density
import ..gravitational_potential
import ..max_time_step
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

function drift!(profile::SphericalProfile, m, h::Real, explap)
    psi = profile.psi

    for field in 1:length(m)
        psi[:, field] .= (explap[:, :, field] * (profile.r .* psi[:, field])) ./ profile.r # u
        # psi[:, field] .= explap[:, :, field] * psi[:, field] # psi (not u)
    end

    return profile
end

function compute_explaps_imag(profile::SphericalProfile, m, h)
    resol = size(profile.psi, 1)
    nfields = size(profile.psi, 2)

    explaps = zeros(resol, resol, nfields)
    for field in 1:nfields
        explaps[:, :, field] .= exp(+h / 2 / m[1, field] * laplacian_u(profile))
    end

    return explaps
end

function max_time_step(profile::SphericalProfile)

    dr = dr_element(profile)

    return dr^2/π
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

function d1_dr1(profile::SphericalProfile)
    resol_r = size(profile.r, 1)

    return d1_dr1(resol_r)
end

function d2_dr2(resol)
    out = Tridiagonal(ones(resol - 1), -2 * ones(resol), ones(resol - 1))
    # Assume derivand has f[begin-1] = f[begin]
    out[begin, begin] = -1.0

    # Asymptote at R=R_max
    out[resol, resol] = -1.0

    # out = Matrix(out)
    # out[begin, 1] = 1
    # out[begin, 2] = -2
    # out[begin, 3] = 1

    # out[end, end] = 1
    # out[end, end-1] = -2
    # out[end, end-2] = 1

    return out
end

function d2_dr2(profile::SphericalProfile)
    resol_r = size(profile.r, 1)

    return d2_dr2(resol_r)
end

function laplacian(profile::SphericalProfile)
    dr = dr_element(profile)

    resol_r = size(profile.r, 1)

    d1 = d1_dr1(profile) / dr
    d1 .*= 2 ./ profile.r
    @assert d1 isa Tridiagonal
    d2 = d2_dr2(profile) / dr^2
    @assert d2 isa Tridiagonal

    return d1 + d2
end

function laplacian_u(profile::SphericalProfile)
    dr = dr_element(profile)

    out = d2_dr2_vanish_r0(profile) / dr^2

    return out
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
