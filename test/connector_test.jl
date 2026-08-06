@testsnippet ConnectorSetup begin
    using AtmosphericDeposition
    using Test, ModelingToolkit, Dates, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using EarthSciData, Aerosol

    domain = DomainInfo(
        DateTime(2016, 2, 1),
        DateTime(2016, 2, 2);
        latrange = deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
        lonrange = deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
        levrange = 1:10
    )
end

# Scalar-landuse `DryDepositionGas` test disabled: that constructor no longer
# exists in `src/dry_deposition.jl` (the exported `DryDepositionGas` now
# resolves to the unrelated Seinfeld-Pandis ch. 19 component). The
# chemistry-side coupling it exercised now runs through
# `DryDepositionGasFractional`, which uses the same `DryDepositionGasCoupler`
# (see the fractional test further below).
#=
@testitem "GasChemExt SuperFast DryDeposition" begin
    using AtmosphericDeposition, GasChem, EarthSciMLBase, ModelingToolkit
    using Test

    model = couple(SuperFast(), DryDepositionGas())
    sys = convert(System, model)
    eqs = string(equations(sys))

    # Verify that GasChem species are coupled to dry deposition rate constants
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_HNO3")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_NO2")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_O3")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_H2O2")
    @test contains(eqs, "SuperFast₊DryDepositionGas_k_HCHO")
end
=#

@testitem "GasChemExt SuperFast WetDeposition" begin
    using AtmosphericDeposition, GasChem, EarthSciMLBase, ModelingToolkit
    using Test

    model = couple(SuperFast(), WetDeposition())
    sys = convert(System, model)
    eqs = string(equations(sys))

    # Verify that GasChem species are coupled to wet deposition rate constants
    @test contains(eqs, "SuperFast₊WetDeposition_k_othergas")
end

@testitem "AerosolExt" setup = [ConnectorSetup] begin
    model = couple(
        GEOSFP("4x5", domain),
        WetDeposition(),
        ElementalCarbon(),
        DryDepositionAerosol()
    )
    sys = convert(System, model)

    eqs = equations(sys)
    @test contains(string(eqs), "ElementalCarbon₊DryDepositionAerosol_k")
    @test contains(string(eqs), "ElementalCarbon₊WetDeposition_k_particle")
end

@testitem "EarthSciDataExt" setup = [ConnectorSetup] begin
    model = couple(
        GEOSFP("4x5", domain),
        WetDeposition(),
        DryDepositionAerosol(),
        ElementalCarbon(),
    )

    sys = convert(System, model)

    eqs = string(observed(sys))
    wanteq = "GEOSFP₊A1₊USTAR(t)"
    @test contains(eqs, wanteq)
    wanteq = "WetDeposition₊cloudFrac(t) ~ GEOSFP₊A3cld₊CLOUD(t)"
    @test contains(eqs, wanteq)
end

@testitem "EarthSciDataExt fractional gas dry deposition" setup = [ConnectorSetup] begin
    # Verify the new couple2 method binds the fractional gas coupler's
    # lon/lat/season to GEOS-FP and that the 9 surface-meteorology
    # bindings carry over from the scalar-landuse coupler.
    model = couple(
        GEOSFP("4x5", domain),
        DryDepositionGasFractional(),
    )
    sys = convert(System, model)
    eqs = string(equations(sys)) * "\n" * string(observed(sys))

    # New (fractional-specific) bindings.
    @test contains(eqs, "DryDepositionGasFractional₊lon")
    @test contains(eqs, "DryDepositionGasFractional₊lat")
    @test contains(eqs, "season_at")
    @test contains(eqs, "GEOSFP₊t_ref")

    # Carried over from the scalar-landuse coupler.
    @test contains(eqs, "GEOSFP₊A1₊USTAR")
    @test contains(eqs, "GEOSFP₊A1₊TS")
    @test contains(eqs, "GEOSFP₊A1₊SWGDN")

    # No scalar landuse parameter survives (the fractional system never
    # declared one — sanity guard against silent regression).
    @test !contains(eqs, "DryDepositionGasFractional₊landuse")
end
