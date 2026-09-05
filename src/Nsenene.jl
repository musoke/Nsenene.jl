module Nsenene

export density
export total_mass, total_masses
export gravitational_potential

export CylindricalProfile
export SphericalProfile

"""
    ground_state(m, Lambda, M, profile_type, resol; r_max=1.0)

Calculate an approximate ground state for particle masses `m`, interactions `Lambda` and total masses `M`.
"""
function ground_state(m, Lambda, M, profile_type, resol; r_max=1.0)
    nfields = size(m)[end]

    profile = profile_type(resol, r_max, nfields)
    r = Nsenene.radius(profile)

    # Gaussian initial condition in each field
    selectdim(profile.psi, ndims(m), 1) .= 100 * exp.(-10r .^ 2)
    selectdim(profile.psi, ndims(m), 2) .= 100 * exp.(-10r .^ 2)

    # Random perturbations
    profile.psi[:, :, :] += randn(size(profile.psi)) ./ r

    relax!(profile, m, Lambda, M)

    return profile
end

function relax!(profile, m, Lambda, M; maxsteps=1000, nsteps=10, atol=0, rtol=1e-3)
    h = max_time_step(profile) * 10

    Nsenene.normalize_mass!(profile, m, M)

    rho_max = Float64[]
    push!(rho_max, maximum(density(profile, m)))

    converged = false

    for step in nsteps:nsteps:maxsteps
        kick_drift_kick!(profile, h, m, Lambda, M, nsteps)

        push!(rho_max, maximum(density(profile, m)))

        if isapprox(rho_max[end], rho_max[end - 1]; atol=atol, rtol=rtol)
            @info "Converged after $step steps"
            converged = true
            break
        end
    end

    if converged
    else
        @error "Relaxation did not converge after $maxsteps steps" rho_max
    end

    return profile
end

function normalize_mass!(profile, m, target_masses)
    current_masses = total_masses(profile, m)

    # Careful of fields with no mass
    mass_ratios = target_masses ./ (current_masses + (current_masses .== 0))

    profile.psi .*= reshape(mass_ratios, size(m)) .^ 0.5

    return profile
end

function kick_drift_kick!(profile, h, m, Lambda, target_masses, nsteps)
    explaps = compute_explaps_imag(profile, m, h)

    for step in 1:nsteps
        kick!(profile, m, Lambda, h / 2)
        drift!(profile, m, h, explaps)
        kick!(profile, m, Lambda, h / 2)

        normalize_mass!(profile, m, target_masses)

        @assert !any(isnan.(profile.psi))
    end
end

"""
    kick(profile, h, m)
"""
function kick!(profile, m, Lambda, h)
    psi = profile.psi

    V = m .* gravitational_potential(profile, m)
    interaction_potential!(V, profile, m, Lambda)
    @assert size(V) == size(profile.psi)

    psi .*= exp.(-h .* V)

    return profile
end

"""
    interaction_potential!(V, profile, m, Lambda)

Add the interaction potential to the potential array V
"""
function interaction_potential! end

function drift! end

function compute_explaps_imag end

function max_time_step end

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

function radius end

include("cylindrical.jl")
include("spherical.jl")

import .Cylindrical: CylindricalProfile
import .Spherical: SphericalProfile

end
