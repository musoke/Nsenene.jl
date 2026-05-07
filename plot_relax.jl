using CairoMakie
using LinearAlgebra
using NPZ
using Nsenene
using ProgressMeter

function add_pyul_soliton!(psi, r, mass)
    profile::Vector{Float64} = npzread(
        joinpath(@__DIR__, "..", "UltraDark", "examples", "initial_f.npy")
    )
    delta_x = 0.00001
    alpha = (mass / 3.883)^2
    beta = 2.454

    for i in eachindex(psi)
        distfromcentre = r[i]
        pyul_idx = trunc(Int, alpha^0.5 * (distfromcentre / delta_x + 1))

        if pyul_idx < length(profile)
            f = alpha * profile[pyul_idx]
            psi[i] += f
        else
            psi[i] += 0
        end
    end
end

function plot_relax()
    resol = 2^9
    r_max = 1e+1
    nfields = 2
    p = SphericalProfile(resol, r_max, nfields)
    m = ones(1, nfields)

    target_masses = zeros(nfields)
    target_masses[1] = 50.0

    add_pyul_soliton!(@view(p.psi[:, 1]), p.r, target_masses[1])
    @show target_masses
    @show total_masses(p, m)

    fig = Figure()
    ax_rho = Axis(fig[1, 1]; xlabel=L"$r$", ylabel=L"$\rho(r)$", xscale=log10, yscale=log10)
    xlims!(ax_rho, 0.9 * p.r[1], r_max)
    ylims!(ax_rho, 1e-5, 1e+5)
    vlines!(ax_rho, p.r[1]; alpha=0.1, color=:black, linestyle=:dot)
    scatter!(ax_rho, [p.r[1]], [abs2(p.psi[1, 1])]; alpha=0.1, color=:black)

    lines!(ax_rho, p.r, abs2.(p.psi[:, 1]); color=:black)

    nits = 10
    nsteps = 1000

    M = Float64[total_mass(p, m)]
    rho_max = Float64[maximum(abs2.(p.psi))]
    @showprogress for i in 1:nits
        phi = gravitational_potential(p, m)
        h = Nsenene.max_time_step(p, phi) / 1000

        Nsenene.kick_drift_kick!(p, h, m, target_masses, nsteps)

        lines!(
            ax_rho,
            p.r,
            abs2.(p.psi[:, 1]);
            linestyle=:dash,
            alpha=0.6,
            color=i,
            colorrange=(0, nits),
        )

        push!(M, total_mass(p, m))
        push!(rho_max, maximum(abs2.(p.psi)))
    end

    ax_rho_max = Axis(
        fig[2, 1];
        xscale=identity,
        yscale=log10,
        xlabel="iteration",
        ylabel=L"$\max(\rho(r))$",
    )
    scatter!(ax_rho_max, rho_max)

    return save("f.pdf", fig)
end

plot_relax()
