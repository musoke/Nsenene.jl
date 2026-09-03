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

function radius(p::SphericalProfile)
    return p.r
end

function radius(p::CylindricalProfile)
    return @. (p.R^2 + p.z^2)^0.5
end

function plot_relax()
    resol = 2^8
    r_max = 1e+0
    nfields = 2
    p = CylindricalProfile(resol, r_max, nfields)
    ps = SphericalProfile(resol, r_max, nfields)

    m = ones(1, 1, nfields)
    ms = ones(1, nfields)
    m[1] = 1
    m[2] = 1
    @show m

    Lambda = zeros(nfields, nfields)
    Lambda[2, 1] = Lambda[1, 2] = 0.1

    target_masses = zeros(nfields)
    target_masses[1] = 90.0
    target_masses[2] = 10.0
    @show target_masses

    r = radius(p)

    fig = Figure()
    ax_rho_1_r = Axis(
        fig[1, 1];
        xlabel=L"$R$",
        ylabel=L"$\rho_1(z=0, R)$",
        xscale=identity,
        yscale=Makie.pseudolog10,
    )
    ax_rho_2_r = Axis(
        fig[2, 1];
        xlabel=L"$R$",
        ylabel=L"$\rho_2(z=0, R)$",
        xscale=identity,
        yscale=Makie.pseudolog10,
    )
    ax_rho_1_z = Axis(
        fig[1, 2];
        xlabel=L"$z$",
        ylabel=L"$\rho_1(z, R=0)$",
        xscale=identity,
        yscale=Makie.pseudolog10,
    )
    ax_rho_2_z = Axis(
        fig[2, 2];
        xlabel=L"$z$",
        ylabel=L"$\rho_2(z, R=0)$",
        xscale=identity,
        yscale=Makie.pseudolog10,
    )

    # xlims!(ax_rho, 0.9 * p.r[1], r_max)
    # ylims!(ax_rho, 1e-0, 1e+5)
    ylims!(ax_rho_1_r, -1, 1e+5)
    ylims!(ax_rho_2_r, -1, 1e+5)
    ylims!(ax_rho_1_z, -1, 1e+5)
    ylims!(ax_rho_2_z, -1, 1e+5)

    # Plot target profile
    if target_masses[1] != 0
        add_pyul_soliton!(@view(p.psi[:, :, 1]), r, target_masses[1])
        lines!(
            ax_rho_1_r,
            p.R[1, :],
            abs2.(p.psi[resol ÷ 2, :, 1]);
            color=:black,
            linewidth=4,
            label="PyUltraLight profile",
        )
        lines!(
            ax_rho_1_z,
            p.z[:, 1],
            abs2.(p.psi[:, 1, 1]);
            color=:black,
            linewidth=4,
            label="PyUltraLight profile",
        )
    end
    if target_masses[2] != 0
        add_pyul_soliton!(@view(p.psi[:, :, 2]), r, target_masses[2])
        lines!(
            ax_rho_2_r,
            p.R[1, :],
            abs2.(p.psi[resol ÷ 2, :, 2]);
            color=:black,
            linewidth=4,
            label="PyUltraLight profile",
        )
        lines!(
            ax_rho_2_z,
            p.z[:, 1],
            abs2.(p.psi[:, 1, 2]);
            color=:black,
            linewidth=4,
            label="PyUltraLight profile",
        )
    end

    # Gaussian initial condition in each field
    selectdim(p.psi, ndims(m), 1) .= 100 * exp.(-10r .^ 2)
    selectdim(p.psi, ndims(m), 2) .= 100 * exp.(-10r .^ 2)
    p.psi[:, :, :] += randn(size(p.psi)) ./ r / 1

    selectdim(ps.psi, ndims(ms), 1) .= 100 * exp.(-10ps.r .^ 2)
    selectdim(ps.psi, ndims(ms), 2) .= 100 * exp.(-10ps.r .^ 2)
    ps.psi[:, :] += randn(size(ps.psi)) ./ ps.r / 1

    # Careful of fields with no mass
    new_masses = total_masses(p, m)
    @show mass_ratios = target_masses ./ (new_masses + (new_masses .== 0))
    p.psi .*= reshape(mass_ratios, size(m)) .^ 0.5

    new_masses = total_masses(ps, ms)
    @show mass_ratios = target_masses ./ (new_masses + (new_masses .== 0))
    ps.psi .*= reshape(mass_ratios, size(ms)) .^ 0.5

    # Plot initial condition
    lines!(
        ax_rho_1_r,
        p.R[1, :],
        abs2.(p.psi[resol ÷ 2, :, 1]);
        color=:grey,
        linewidth=4,
        label="Initial profile",
    )
    lines!(
        ax_rho_1_z,
        p.z[:, 1],
        abs2.(p.psi[:, 1, 1]);
        color=:grey,
        linewidth=4,
        label="Initial profile",
    )

    lines!(
        ax_rho_2_r,
        p.R[1, :],
        abs2.(p.psi[resol ÷ 2, :, 2]);
        color=:grey,
        linewidth=4,
        label="Initial profile",
    )
    lines!(
        ax_rho_2_z,
        p.z[:, 1],
        abs2.(p.psi[:, 1, 2]);
        color=:grey,
        linewidth=4,
        label="Initial profile",
    )

    @show target_masses
    @show total_masses(p, m)

    nits = 50
    nsteps = 10

    M = Float64[total_mass(p, m)]
    rho_max = Float64[maximum(abs2.(p.psi))]
    @show h = Nsenene.max_time_step(p) * 100
    @show hs = Nsenene.max_time_step(ps) * 100

    @showprogress for i in 1:nits
        Nsenene.kick_drift_kick!(p, h, m, Lambda, target_masses, nsteps)
        Nsenene.kick_drift_kick!(ps, hs, ms, Lambda, target_masses, nsteps)

        lines!(
            ax_rho_1_r,
            p.R[1, :],
            abs2.(p.psi[resol ÷ 2, :, 1]);
            linestyle=:dash,
            alpha=0.6,
            color=i,
            colorrange=(0, nits),
        )

        lines!(
            ax_rho_1_z,
            p.z[:, 1],
            abs2.(p.psi[:, 1, 1]);
            linestyle=:dash,
            alpha=0.6,
            color=i,
            colorrange=(0, nits),
        )

        lines!(
            ax_rho_2_r,
            p.R[1, :],
            abs2.(p.psi[resol ÷ 2, :, 2]);
            linestyle=:dash,
            alpha=0.6,
            color=i,
            colorrange=(0, nits),
        )

        lines!(
            ax_rho_2_z,
            p.z[:, 1],
            abs2.(p.psi[:, 1, 2]);
            linestyle=:dash,
            alpha=0.6,
            color=i,
            colorrange=(0, nits),
        )

        push!(M, total_mass(p, m))
        push!(rho_max, maximum(abs2.(p.psi)))
    end

    ax_rho_1 = Axis(fig[3, 1]; xlabel=L"z", ylabel=L"R", aspect=DataAspect())
    ax_rho_2 = Axis(fig[3, 2]; xlabel=L"z", ylabel=L"R", aspect=DataAspect())

    ylims!(ax_rho_1, 0, r_max / 2)
    ylims!(ax_rho_2, 0, r_max / 2)

    heatmap!(ax_rho_1, p.z[:, 1], p.R[1, :], abs2.(p.psi[:, :, 1]))
    heatmap!(ax_rho_2, p.z[:, 1], p.R[1, :], abs2.(p.psi[:, :, 2]))

    # ax_rho_max = Axis(
    #     fig[4, 1:2];
    #     xscale=identity,
    #     yscale=Makie.pseudolog10,
    #     xlabel="iteration",
    #     ylabel=L"$\max(\rho)$",
    # )
    # scatter!(ax_rho_max, rho_max)

    save("cyl-m=$m-target_masses=$target_masses-r_max=$r_max-rep.pdf", fig)
    # save("cyl-m=$m-target_masses=$target_masses-h=$h-r_max=$r_max.pdf", fig)

    fig_comp = Figure(; title=L"$m=%$m$, $\Lambda = %$Lambda$")
    ax1 = Axis(fig_comp[1, 1]; yscale=Makie.pseudolog10)
    ax2 = Axis(fig_comp[2, 1]; yscale=Makie.pseudolog10)

    lines!(ax1, ps.r, abs2.(ps.psi[:, 1]); label="spherical")
    lines!(
        ax1,
        p.R[1, :],
        abs2.(p.psi[resol ÷ 2, :, 1]);
        linestyle=:dash,
        label="cylindrical, R",
    )
    lines!(ax1, p.z[:, 1], abs2.(p.psi[:, 1, 1]); linestyle=:dot, label="cylindrical, z")

    lines!(ax2, ps.r, abs2.(ps.psi[:, 2]); label="spherical")
    lines!(
        ax2,
        p.R[1, :],
        abs2.(p.psi[resol ÷ 2, :, 2]);
        linestyle=:dash,
        label="cylindrical, R",
    )
    lines!(ax2, p.z[:, 1], abs2.(p.psi[:, 1, 2]); linestyle=:dot, label="cylindrical, z")

    fig_comp[1:2, 2] = Legend(fig_comp, ax1, ""; framevisible=false)

    save("compare.pdf", fig_comp)
    # save("compare-m=$m-target_masses=$target_masses-r_max=$r_max-rep.pdf", fig_comp)

    return p
end

plot_relax()
