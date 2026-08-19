@testsnippet FractionalSetup begin
    using AtmosphericDeposition
    using AtmosphericDeposition: FractionalWesleyVd, landuse_frac_at,
        DryDepGasFractional, reload_landuse_fractions!, set_landuse_data_path!
    using Test, DynamicQuantities, ModelingToolkit
    using StaticArrays
    using NCDatasets
end

@testitem "landuse_frac_at bundled CONUS fractions" setup=[FractionalSetup] begin
    # Bundled NetCDF is MODIS MCD12C1 v6.1 re-classed to Wesely classes.
    # Verify the physical sanity of fractions at well-known locations:
    # every cell must sum to 1.0, and the dominant class must match the
    # geography. Tolerances are loose — we want catastrophic-bug regression
    # protection, not a pixel-exact comparison to MODIS.

    function frac(lon_deg, lat_deg)
        [landuse_frac_at(deg2rad(lon_deg), deg2rad(lat_deg), i) for i in 1:11]
    end

    # All cells in the CONUS interior sum to exactly 1.0 (the area weighting in
    # FractionalWesleyVd relies on this; any leftover is silently dropped).
    for (lon, lat) in [(-100.0, 40.0), (-90.0, 35.0), (-118.0, 34.0), (-75.0, 40.0)]
        @test sum(frac(lon, lat)) ≈ 1.0 atol = 1.0e-5
    end

    # Class-1 (urban) dominant over LA basin.
    @test frac(-118.25, 34.05)[1] > 0.5
    # Class-2 (agricultural) dominant over Iowa.
    @test frac(-93.5, 42.0)[2] > 0.5
    # Class-7 (water) ≈ 1.0 over Gulf of Mexico.
    @test frac(-88.0, 27.0)[7] > 0.95
end

@testitem "landuse_frac_at lookup varies with lat/lon" setup=[FractionalSetup] begin
    # Synthesise a NetCDF with class 1 = 1.0 in the western half, class 7 =
    # 1.0 in the eastern half, so the lookup must actually use lon. Bin to
    # the same grid as the bundled file (96 × 51 cells).
    mktempdir() do dir
        path = joinpath(dir, "split.nc")
        n_lon, n_lat, n_class = 96, 51, 11
        lon = collect(Float32, range(-125.625f0; step = 0.625f0, length = n_lon))
        lat = collect(Float32, range(24.5f0; step = 0.5f0, length = n_lat))
        frac = zeros(Float32, n_lon, n_lat, n_class)
        for j in 1:n_lat, i in 1:n_lon
            if lon[i] < -96.0
                frac[i, j, 1] = 1.0   # urban (western half)
            else
                frac[i, j, 7] = 1.0   # water (eastern half)
            end
        end
        NCDataset(path, "c") do ds
            defDim(ds, "lon", n_lon); defDim(ds, "lat", n_lat); defDim(ds, "class", n_class)
            defVar(ds, "lon", lon, ("lon",))
            defVar(ds, "lat", lat, ("lat",))
            defVar(ds, "fraction", frac, ("lon", "lat", "class"))
        end

        old_cache = nothing
        try
            set_landuse_data_path!(path)
            # Western point — urban
            @test landuse_frac_at(deg2rad(-120.0), deg2rad(40.0), 1) == 1.0
            @test landuse_frac_at(deg2rad(-120.0), deg2rad(40.0), 7) == 0.0
            # Eastern point — water
            @test landuse_frac_at(deg2rad(-70.0), deg2rad(40.0), 1) == 0.0
            @test landuse_frac_at(deg2rad(-70.0), deg2rad(40.0), 7) == 1.0
            # Out-of-range query (south of grid) clamps to nearest cell.
            @test landuse_frac_at(deg2rad(-120.0), deg2rad(10.0), 1) == 1.0
        finally
            # Reset to the bundled default for subsequent tests.
            default_path = joinpath(@__DIR__, "..", "data", "landuse_fractions_conus.nc")
            set_landuse_data_path!(default_path)
        end
    end
end

@testitem "FractionalWesleyVd single-class equality" setup=[FractionalSetup] begin
    # For each Wesely class i, putting all weight on i must give exactly the
    # single-class deposition velocity 1/(RaRb + Rc_i). This is the degenerate
    # case in which the area-weighting cancels entirely.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    RaRb = 50.0

    for i in 1:11
        fractions = zeros(11)
        fractions[i] = 1.0

        rc_single = WesleySurfaceResistance(
            AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, i,
            false, false, false, true
        )
        vd_frac = FractionalWesleyVd(
            RaRb, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
            fractions...,
            false, false, false, true
        )
        @test vd_frac ≈ 1 / (RaRb + rc_single)
    end
end

@testitem "FractionalWesleyVd two-class arithmetic" setup=[FractionalSetup] begin
    # 50% agricultural, 50% mixed forest — area-weighted Vd check for O3.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    RaRb = 50.0
    rc_ag = WesleySurfaceResistance(
        AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, 2,
        false, false, false, true
    )
    rc_mf = WesleySurfaceResistance(
        AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, 6,
        false, false, false, true
    )
    vd_expected = 0.5 / (RaRb + rc_ag) + 0.5 / (RaRb + rc_mf)

    vd_frac = FractionalWesleyVd(
        RaRb, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
        0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0,
        false, false, false, true
    )
    @test vd_frac ≈ vd_expected
end

@testitem "FractionalWesleyVd three-class arithmetic" setup=[FractionalSetup] begin
    # Realistic mixed cell for HNO3: 50% deciduous, 30% coniferous, 20% urban.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    RaRb = 50.0
    rc_dec = WesleySurfaceResistance(
        AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason, 4,
        false, false, false, false
    )
    rc_con = WesleySurfaceResistance(
        AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason, 5,
        false, false, false, false
    )
    rc_urb = WesleySurfaceResistance(
        AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason, 1,
        false, false, false, false
    )
    vd_expected = 0.5 / (RaRb + rc_dec) + 0.3 / (RaRb + rc_con) +
        0.2 / (RaRb + rc_urb)

    vd_frac = FractionalWesleyVd(
        RaRb, AtmosphericDeposition.HNO3Data, G, Ts, θ, iSeason,
        0.2, 0.0, 0.0, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        false, false, false, false
    )
    @test vd_frac ≈ vd_expected
end

@testitem "FractionalWesleyVd is below the Rc-blending form" setup=[FractionalSetup] begin
    # Direction guard. Blending 1/Rc first and adding a single Ra+Rb afterwards
    # is the Ra+Rb -> 0 limit of the correct expression, and by concavity of
    # c -> c/(1+ac) it always gives Vd >= this one, with equality only when
    # every class shares one Rc. A regression to that form must fail here.
    G, Ts, θ, iSeason = 800.0, 25.0, 0.0, 1
    # 96% water + 4% wetland: the widest Rc spread in the bundled CONUS data.
    fr = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.96, 0.0, 0.04, 0.0, 0.0)
    rc = [WesleySurfaceResistance(AtmosphericDeposition.O3Data, G, Ts, θ, iSeason, i,
        false, false, false, true) for i in 1:11]

    for RaRb in (10.0, 50.0, 100.0, 300.0, 500.0)
        vd_new = FractionalWesleyVd(
            RaRb, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
            fr..., false, false, false, true
        )
        vd_old = 1 / (RaRb + 1 / sum(fr[i] / rc[i] for i in 1:11))
        @test vd_new < vd_old
    end

    # Equal Rc across the populated classes => the two forms coincide.
    fr_eq = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0)  # water, barren
    if rc[7] ≈ rc[8]
        vd_new = FractionalWesleyVd(
            100.0, AtmosphericDeposition.O3Data, G, Ts, θ, iSeason,
            fr_eq..., false, false, false, true
        )
        @test vd_new ≈ 1 / (100.0 + rc[7])
    end
end

@testitem "DryDepGasFractional unit" setup=[FractionalSetup] begin
    @parameters T [unit = u"K"]
    @parameters z [unit = u"m"]
    @parameters z₀ [unit = u"m"]
    @parameters u_star [unit = u"m/s"]
    @parameters L [unit = u"m"]
    @parameters ρA [unit = u"kg*m^-3"]
    @parameters G [unit = u"W*m^-2"]
    @parameters θ
    @parameters iSeason, lev

    # Concrete fractions: f_mixedforest = 1.0, others = 0.0.
    # GasData overload — this is the active path used by the
    # `DryDepositionGasFractional` broadcast (via `datas = [...]`).
    @test ModelingToolkit.get_unit(
        DryDepGasFractional(
            lev, z, z₀, u_star, L, ρA,
            AtmosphericDeposition.So2Data,
            G, T, θ, iSeason,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            false, false, true, false
        )
    ) == u"m/s"
end

@testitem "DryDepositionGasFractional constructor" setup=[FractionalSetup] begin
    sys = DryDepositionGasFractional()

    @test sys isa ModelingToolkit.AbstractSystem

    # Fractions are now ordinary parameters (set per cell by the coupler),
    # not system unknowns. This matches the scalar-landuse baseline's
    # symbolic shape: only v_* / k_* are unknowns; everything else is a
    # parameter. lon/lat live on the GEOS-FP side and aren't dep-system
    # parameters anymore.
    param_str = string.(parameters(sys))
    @test any(contains.(param_str, "f_urban"))
    @test any(contains.(param_str, "f_agricultural"))
    @test any(contains.(param_str, "f_mixedforest"))
    @test any(contains.(param_str, "f_rockyshrubs"))
    @test any(contains.(param_str, "season"))
    @test !any(contains.(param_str, "landuse"))
    @test !any(contains.(param_str, "lon"))
    @test !any(contains.(param_str, "lat"))

    var_str = string.(unknowns(sys))
    @test !any(contains.(var_str, "f_urban"))
end
