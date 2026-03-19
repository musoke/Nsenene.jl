module spherical

struct SphericalProfile
    r::Matrix{Float16}
    psi::Matrix{Complex{Float16}}
end

function SphericalProfile(resol::Integer, length::Real, n::Integer)
    r = range(length / (resol + 1), length, resol)  # FIXME: is this the best way to address r=0?
    # For easier broadcasting later
    r = reshape(r, (1, resol))

    psi = zeros(Complex{Float16}, n, resol)

    return SphericalProfile(r, psi)
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
    return sum(densities(profile, m); dims=1)
end

"""
    total_masses(r, psi, m)

Calculate the mass contained by each field with radial coordinates `r` and particle masses `m`.

# Examples

```jldoctest
import Nsenene.spherical: total_masses

N = 1000

m = [1.0, 1.0]
r = range(0, 10, N)
psi = zeros(2, N)

# Ball of radius 1
psi[1, r.<1] .= 1

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

    r = reshape(r, (1, resol))

    @assert size(psi, 2) === size(r, 2)
    @assert size(psi, 1) === size(m, 1)

    dr = diff(r; dims=2)[1]
    M = sum(4π * r .^ 2 .* abs2.(psi); dims=2) .* m * dr

    # @assert size(m) === size(M)

    return M[:, 1]
end

"""
    total_masses(profile::SphericalProfile, m)

Calculate the mass contained by each field with radial coordinates `r` and particle masses `m`.

# Examples

```jldoctest
import Nsenene.spherical: SphericalProfile, total_masses

n = 3
m = ones(n)

p = SphericalProfile(1000, 3.0, n)

# Ball of radius 1
p.psi[1, :] .= p.r[1, :] .< 1
p.psi[2, :] .= p.r[1, :] .< 2

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

end
