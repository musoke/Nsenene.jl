import ..Cylindrical: d1_dR1, dR_element
import ..Cylindrical: d1_dz1, dz_element
import ..Spherical: d1_dr1, dr_element

function energy_density_gravity(profile, m)
    Phi = gravitational_potential(profile, m)
    rho = density(profile, m)

    return 0.5 * Phi .* rho
end

function energy_density_kq(profile::CylindricalProfile, m)
    nfields = length(m)

    grad_R = d1_dR1(profile) ./ dR_element(profile)
    grad_z = d1_dz1(profile) ./ dz_element(profile)

    psi = profile.psi

    out = zeros(Float64, size(profile.psi)[1:2])

    for field in 1:nfields
        out += 0.5 * abs2.(grad_z * psi[:, :, field]) / m[field]
        out += 0.5 * abs2.(psi[:, :, field] * transpose(grad_R)) / m[field]
    end

    return out
end

function energy_density_kq(profile::SphericalProfile, m)
    grad = d1_dr1(profile) ./ dr_element(profile)

    e_kq_per_field = 0.5 * abs2.(grad * profile.psi) ./ m
    e_kq = sum(e_kq_per_field; dims=2)

    return dropdims(e_kq; dims=2)
end

function energy_density_interaction(profile::CylindricalProfile, m, Lambda)
    out = zeros(Float64, size(profile.psi)[1:2])
    nfields = length(m)

    for field_1 in 1:nfields
        for field_2 in field_1:nfields
            psi1 = profile.psi[:, :, field_1]
            psi2 = profile.psi[:, :, field_2]
            out +=
                Lambda[field_1, field_2] / 2 / m[field_1] / m[field_2] * abs2.(psi1) .*
                abs2.(psi2)
        end
    end

    return out
end

function energy_density_interaction(profile::SphericalProfile, m, Lambda)
    out = zeros(Float64, size(profile.psi)[1:1])
    nfields = length(m)

    for field_1 in 1:nfields
        for field_2 in field_1:nfields
            psi1 = profile.psi[:, field_1]
            psi2 = profile.psi[:, field_2]
            out +=
                Lambda[field_1, field_2] / 2 / m[field_1] / m[field_2] * abs2.(psi1) .*
                abs2.(psi2)
        end
    end

    return out
end

function energy_density(profile, m, Lambda)
    out = energy_density_gravity(profile, m)
    out += energy_density_kq(profile, m)
    out += energy_density_interaction(profile, m, Lambda)

    return out
end

function energy(profile, m, Lambda)
    e_density = energy_density(profile, m, Lambda)

    return integrate(e_density, profile)
end

function integrate(integrand, profile::CylindricalProfile)
    R = profile.R
    dR = dR_element(profile)
    dz = dz_element(profile)

    return sum(2π * R .* integrand) * dR * dz
end

function integrate(integrand, profile::SphericalProfile)
    r = profile.r
    dr = dr_element(profile)

    return sum(4π * r .^ 2 .* integrand) * dr
end
