using NCDatasets
using Dates: unix2datetime, month

export reload_landuse_fractions!, set_landuse_data_path!

# Default path to the bundled fractions NetCDF.
const _DEFAULT_LANDUSE_NC_PATH = joinpath(@__DIR__, "..", "data", "landuse_fractions_conus.nc")
const _LANDUSE_NC_PATH = Ref{String}(_DEFAULT_LANDUSE_NC_PATH)

# Concretely-typed cache fields, populated eagerly by `__init__()`. This mirrors
# the GasChem FastJX pattern (e.g. `const top_flux::SVector{18, Float32}`) so
# `landuse_frac_at` is type-stable and allocation-free in the hot RHS path.
# A `Ref{Union{Nothing, NamedTuple{...}}}` would force every cache access to
# narrow the union at runtime, causing per-call boxing and GC pressure
# (~2.8× slowdown observed in 9-day grid runs prior to this change).
const _LANDUSE_FRAC     = Ref{Array{Float32, 3}}()  # (lon, lat, class)
const _LANDUSE_LON_MIN  = Ref{Float64}()
const _LANDUSE_LON_STEP = Ref{Float64}()
const _LANDUSE_LAT_MIN  = Ref{Float64}()
const _LANDUSE_LAT_STEP = Ref{Float64}()
const _LANDUSE_N_LON    = Ref{Int}()
const _LANDUSE_N_LAT    = Ref{Int}()

function _populate_landuse_cache!()
    path = _LANDUSE_NC_PATH[]
    isfile(path) || error(
        "Land-use fractions NetCDF not found at $path. " *
            "Generate it with `julia --project=data data/make_placeholder.jl` " *
            "(placeholder) or `julia --project=data data/preprocess_landuse.jl` " *
            "(real MODIS-derived fractions; requires NASA Earthdata credentials).",
    )
    ds = NCDataset(path)
    lon = Float32.(ds["lon"][:])
    lat = Float32.(ds["lat"][:])
    frac = Float32.(ds["fraction"][:, :, :])
    close(ds)

    _LANDUSE_FRAC[]     = frac
    _LANDUSE_LON_MIN[]  = Float64(lon[1])
    _LANDUSE_LON_STEP[] = Float64(lon[2]) - Float64(lon[1])
    _LANDUSE_LAT_MIN[]  = Float64(lat[1])
    _LANDUSE_LAT_STEP[] = Float64(lat[2]) - Float64(lat[1])
    _LANDUSE_N_LON[]    = length(lon)
    _LANDUSE_N_LAT[]    = length(lat)
    return nothing
end

# Eager load at package init so the first RHS call doesn't trigger I/O or
# first-call branching. If the NetCDF is missing or malformed, warn and leave
# the refs unassigned: `landuse_frac_at` will then surface an `UndefRefError`
# at first call, paired with the warning here for diagnostic context.
function __init__()
    try
        _populate_landuse_cache!()
    catch err
        @warn "Land-use fractions NetCDF could not be loaded; \
               `landuse_frac_at` will error on first use." err
    end
end

"""
    reload_landuse_fractions!()

Re-read the NetCDF at the current path into the cache. Useful after
regenerating the file with `preprocess_landuse.jl`.
"""
function reload_landuse_fractions!()
    _populate_landuse_cache!()
    return nothing
end

"""
    set_landuse_data_path!(path)

Override the NetCDF path used by `landuse_frac_at` and reload immediately.
Primarily used by tests to point at a synthetic fixture.
"""
function set_landuse_data_path!(path::AbstractString)
    _LANDUSE_NC_PATH[] = String(path)
    _populate_landuse_cache!()
    return nothing
end

# Area fraction of Wesely (1989) land-use class `class_index` at the
# given longitude/latitude (RADIANS, matching `DomainInfo`). Indices
# follow the `wesleyLandUse` enum:
#   1 = urban, 2 = agricultural, 3 = range, 4 = deciduous, 5 = coniferous,
#   6 = mixed forest, 7 = water, 8 = barren, 9 = wetland, 10 = rangeag,
#   11 = rocky shrubs.
#
# Nearest-cell lookup on the uniform CONUS grid bundled in
# `data/landuse_fractions_conus.nc`. Query points outside the grid are
# clamped to the nearest boundary cell.
function landuse_frac_at(lon_rad, lat_rad, class_index)::Float64
    lon_deg = rad2deg(lon_rad)
    lat_deg = rad2deg(lat_rad)
    i_lon = clamp(round(Int, (lon_deg - _LANDUSE_LON_MIN[]) / _LANDUSE_LON_STEP[]) + 1, 1, _LANDUSE_N_LON[])
    i_lat = clamp(round(Int, (lat_deg - _LANDUSE_LAT_MIN[]) / _LANDUSE_LAT_STEP[]) + 1, 1, _LANDUSE_N_LAT[])
    @inbounds Float64(_LANDUSE_FRAC[][i_lon, i_lat, class_index])
end
# Registered as symbolic so that `lon`/`lat` are preserved as parameter
# references in the System equations (otherwise the symbolic engine
# might fold based on the lookup logic).
@register_symbolic landuse_frac_at(lon, lat, class_index)
landuse_frac_at(::DynamicQuantities.Quantity, ::DynamicQuantities.Quantity, ::Any) = 0.0
ModelingToolkit.get_unit(::typeof(landuse_frac_at)) = 1.0

# Wesely (1989) season index from absolute Unix time, assuming
# Northern-Hemisphere mid-latitudes:
#   Jun–Aug → 1 (Midsummer)
#   Sep–Oct → 2 (Autumn)
#   Nov     → 3 (LateAutumn)
#   Dec–Feb → 4 (Winter)
#   Mar–May → 5 (Transitional)
#
# Called from `couple2(::DryDepositionGasCoupler, ::GEOSFPCoupler)` in
# `ext/EarthSciDataExt.jl` as `season_at(gp.t_ref + t)`. For
# Southern-Hemisphere or tropical work, override `season` manually instead of
# binding to this function.
#
# TODO: lat-aware seasons for global use.
function season_at(t_abs_seconds)
    m = month(unix2datetime(t_abs_seconds))
    return m == 6 || m == 7 || m == 8  ? Int(wesleyMidsummer)    :
           m == 9 || m == 10           ? Int(wesleyAutumn)       :
           m == 11                     ? Int(wesleyLateAutumn)   :
           m == 12 || m == 1 || m == 2 ? Int(wesleyWinter)       :
                                         Int(wesleyTransitional)
end
@register_symbolic season_at(t_abs_seconds)
season_at(::DynamicQuantities.Quantity) = 1
ModelingToolkit.get_unit(::typeof(season_at)) = 1.0
