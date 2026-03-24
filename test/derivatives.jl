import Nsenene: Cylindrical, CylindricalProfile, Spherical, SphericalProfile

import LinearAlgebra: Tridiagonal

resol = 32
length = 3.0
nfields = 1
m_ = 1:nfields

@testset "CylindricalProfile" begin
    p = CylindricalProfile(resol, length, nfields)

    dR1 = Cylindrical.d1_dR1(p)
    @test dR1 isa Tridiagonal
    @test size(dR1) == (resol, resol)

    @testset "first derivative of r should be 1" begin
        d = dR1 * p.R / Cylindrical.dR_element(p)

        @test d[begin:end] ≈ ones(resol) atol = 1e-10
    end

    dR2 = Cylindrical.d2_dR2(p)

    @test dR2 isa Tridiagonal
    @test size(dR2) == (resol, resol)

    @testset "second derivative of r should vanish" begin
        d = dR2 * p.R / Cylindrical.dR_element(p)^2

        @test d[(begin + 1):(end - 1)] ≈ zeros(resol - 2) atol = 1e-10
        @test_broken d[begin] ≈ 0.0 atol = 1e-10
        @test_broken d[end] ≈ 0.0 atol = 1e-10
    end

    @testset "second derivative of r^2 should be 2" begin
        d = dR2 * (p.R .^ 2) / Cylindrical.dR_element(p)^2

        @test d[(begin + 1):(end - 1)] ≈ 2 * ones(resol - 2) atol = 1e-10
        @test_broken d[begin] ≈ 2.0 atol = 1e-10
        @test_broken d[end] ≈ 2.0 atol = 1e-10
    end
end

@testset "SphericalProfile" begin
    p = SphericalProfile(resol, length, nfields)

    dr2 = Spherical.d2_dr2(p)

    @test dr2 isa Tridiagonal
    @test size(dr2) == (resol, resol)

    @testset "second derivative of r should vanish" begin
        d = dr2 * p.r / Spherical.dr_element(p)^2

        @test d[(begin + 1):(end - 1)] ≈ zeros(resol - 2) atol = 1e-10
        @test_broken d[begin] ≈ 0.0 atol = 1e-10
        @test_broken d[end] ≈ 0.0 atol = 1e-10
    end

    @testset "second derivative of r^2 should be 2" begin
        d = dr2 * (p.r .^ 2) / Spherical.dr_element(p)^2

        @test d[(begin + 1):(end - 1)] ≈ 2 * ones(resol - 2) atol = 1e-10
        @test_broken d[begin] ≈ 2.0 atol = 1e-10
        @test_broken d[end] ≈ 2.0 atol = 1e-10
    end
end
