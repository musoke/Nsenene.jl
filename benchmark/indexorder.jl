using AbstractFFTs
using BenchmarkTools
using FFTW
using LinearAlgebra

suite = BenchmarkGroup()

suite["spherical"] = BenchmarkGroup()
begin
    resol_r = 10_000
    nfields = 1

    a_long_first = rand(resol_r, nfields)
    r_long_first = range(1.0, 10.0, resol_r)
    suite["spherical"]["r first"] = @benchmarkable abs2.($a_long_first) ./ $r_long_first

    a_long_second = rand(nfields, resol_r)
    r_long_second = reshape(r_long_first, 1, resol_r)
    suite["spherical"]["field first"] = @benchmarkable abs2.($a_long_second) ./
        $r_long_second
end

begin
    resol_R = 1024
    resol_z = 2 * resol_R

    rho_R_first = randn(resol_R, resol_z)
    rho_z_first = permutedims(rho_R_first, (2, 1))

    @assert size(rho_R_first) == (resol_R, resol_z)
    @assert size(rho_z_first) == (resol_z, resol_R)

    fftplan_R_first = AbstractFFTs.plan_rfft(rho_R_first, 2)
    fftplan_z_first = AbstractFFTs.plan_rfft(rho_z_first, 1)

    # fftplan_R_first_patient = AbstractFFTs.plan_rfft(copy(rho_R_first), 2; flags=FFTW.PATIENT)
    # fftplan_z_first_patient = AbstractFFTs.plan_rfft(copy(rho_z_first), 1; flags=FFTW.PATIENT)

    suite["cylindrical"]["FFT"] = BenchmarkGroup()

    suite["cylindrical"]["FFT"]["R first estimate"] = @benchmarkable $fftplan_R_first *
        $rho_R_first
    suite["cylindrical"]["FFT"]["z first estimate"] = @benchmarkable $fftplan_z_first *
        $rho_z_first

    # suite["cylindrical"]["FFT"]["R first patient"] = @benchmarkable $fftplan_R_first_patient * $rho_R_first
    # suite["cylindrical"]["FFT"]["z first patient"] = @benchmarkable $fftplan_z_first_patient * $rho_z_first

    suite["cylindrical"]["FFT"]["R first unplanned"] = @benchmarkable rfft(rho_R_first, 2)
    suite["cylindrical"]["FFT"]["z first unplanned"] = @benchmarkable rfft(rho_z_first, 1)

    rhok_R_first = fftplan_R_first * rho_R_first
    rhok_z_first = fftplan_z_first * rho_z_first

    ifftplan_R_first = AbstractFFTs.plan_irfft(rhok_R_first, resol_z, 2)
    ifftplan_z_first = AbstractFFTs.plan_irfft(rhok_z_first, resol_z, 1)

    suite["cylindrical"]["IFFT"]["R first estimate"] = @benchmarkable $ifftplan_R_first *
        $rhok_R_first
    suite["cylindrical"]["IFFT"]["z first estimate"] = @benchmarkable $ifftplan_z_first *
        $rhok_z_first

    kz = AbstractFFTs.rfftfreq(resol_z)
    D = Tridiagonal(randn(resol_R, resol_R))

    function operator_inverse_R_first(rhok, D, kz, Phik)
        for (i, kz) in enumerate(kz)
            Dk = D - UniformScaling(kz^2)
            Dk[end, end - 1] = 0
            Dk[end, end] = 1

            Phik[:, i] .= Dk \ rhok[:, i]
        end

        return Phik
    end

    function operator_inverse_z_first(rhok, D, kz, Phik)
        for (i, kz) in enumerate(kz)
            Dk = D - UniformScaling(kz^2)
            Dk[end, end - 1] = 0
            Dk[end, end] = 1

            Phik[i, :] .= Dk \ rhok[i, :]
        end

        return Phik
    end

    inv_R_first = operator_inverse_R_first(rhok_R_first, D, kz, similar(rhok_R_first))
    inv_z_first = operator_inverse_z_first(rhok_z_first, D, kz, similar(rhok_z_first))

    @assert all((inv_R_first ≈ transpose(inv_z_first)))

    suite["cylindrical"]["invert tridiagonal"]["R first"] = @benchmarkable operator_inverse_R_first(
        rhok_R_first, D, kz, similar(rhok_R_first)
    )
    suite["cylindrical"]["invert tridiagonal"]["z first"] = @benchmarkable operator_inverse_z_first(
        rhok_z_first, D, kz, similar(rhok_z_first)
    )
end

results = run(suite)

@show judge(
    minimum(results["spherical"]["r first"]), minimum(results["spherical"]["field first"])
)

results
