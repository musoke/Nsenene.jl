module Cylindrical

import LinearAlgebra: UniformScaling, Tridiagonal

using AbstractFFTs: AbstractFFTs
using FFTW: FFTW

import ..densities
import ..density
import ..gravitational_potential
import ..total_masses
import ..total_mass

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
    psi::Array{Complex{Float64},3}
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

function dR_element(profile::CylindricalProfile)
    R = profile.R
    return diff(R; dims=1)[1]
end

function dz_element(profile::CylindricalProfile)
    z = profile.z
    return diff(z; dims=2)[1]
end

function total_masses(profile::CylindricalProfile, m)
    @assert size(profile.psi, 3) === size(m, 3)

    R = profile.R
    z = profile.z
    psi = profile.psi

    dR = dR_element(profile)
    dz = dz_element(profile)

    M = sum(R .* abs2.(psi); dims=(1, 2)) .* m * 2pi * dR * dz

    @assert size(m) === size(M)

    return reshape(M, :)
end

function d1_dR1(resol_R)
    out = 0.5 * Tridiagonal(-ones(resol_R - 1), zeros(resol_R), ones(resol_R - 1))

    # # Asymptote at R=0
    # out[begin, begin] = 0.5

    # # Asymptote at R=end
    # out[end, end] = 0.5

    # Forward difference at R=0
    out[begin, begin] = -1
    out[begin, begin + 1] = 1

    # Backward difference at R=end
    out[end, end - 1] = -1
    out[end, end] = 1

    return out
end

function d1_dR1(profile::CylindricalProfile)
    resol_R = size(profile.R, 1)

    return d1_dR1(resol_R)
end

function d2_dR2(resol_R)
    out = Tridiagonal(ones(resol_R - 1), -2 * ones(resol_R), ones(resol_R - 1))

    # Assymptote at R=0
    out[begin, begin] = -1.0

    # Assymptote at R=R_max
    out[end, end] = -1.0

    return out
end

function d2_dR2(profile::CylindricalProfile)
    resol_R = size(profile.R, 1)

    return d2_dR2(resol_R)
end

function gravitational_potential(profile::CylindricalProfile, m)
    dR = dR_element(profile)
    R = profile.R
    resol_R = size(R, 1)

    dz = dz_element(profile)
    z = profile.z
    resol_z = size(z, 2)
    kz = AbstractFFTs.rfftfreq(resol_z, 2pi / dz)

    M = total_mass(profile, m)
    rho = density(profile, m)

    fftz = AbstractFFTs.plan_rfft(rho, 2)

    rho_Rk = fftz * rho
    rho_Rk .*= 4 * pi * G

    Phi_Rk = similar(rho_Rk)

    phi_boundary = -G * M ./ sqrt.(R[end]^2 .+ z[1, :] .^ 2)
    phi_boundary_k = AbstractFFTs.rfft(phi_boundary, 1)

    _D = Tridiagonal(1 ./ R .* d1_dR1(profile) / dR) + Tridiagonal(d2_dR2(profile) / dR^2)

    for (i, kz) in enumerate(kz)
        D = _D - UniformScaling(kz^2)

        # Required for reasonable speed in the inversion
        # See https://github.com/JuliaLang/LinearAlgebra.jl/issues/1543
        @assert D isa Tridiagonal

        # Enforce boundary condition on Phi_Rk at R=Rmax
        D[end, end - 1] = 0
        D[end, end] = 1
        rho_Rk[end, i] = phi_boundary_k[i]

        Phi_Rk[:, i] .= D \ rho_Rk[:, i]
    end

    Phi = AbstractFFTs.irfft(Phi_Rk, resol_z, 2)

    return Phi
end

end
