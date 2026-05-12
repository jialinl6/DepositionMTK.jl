using NCDatasets
using Dates: unix2datetime, month

export reload_landuse_fractions!, set_landuse_data_path!

# Default path to the bundled fractions NetCDF.
const _DEFAULT_LANDUSE_NC_PATH = joinpath(@__DIR__, "..", "data", "landuse_fractions_conus.nc")
const _LANDUSE_NC_PATH = Ref{String}(_DEFAULT_LANDUSE_NC_PATH)

# Lazily-loaded, cached read of the fractions NetCDF.
const _LANDUSE_CACHE = Ref{Union{Nothing,
    NamedTuple{
        (:lon, :lat, :frac, :lon_min, :lon_step, :lat_min, :lat_step,
            :n_lon, :n_lat, :n_class),
        Tuple{Vector{Float32}, Vector{Float32}, Array{Float32, 3},
            Float64, Float64, Float64, Float64, Int, Int, Int},
    },
}}(nothing)

function _load_landuse()
    cache = _LANDUSE_CACHE[]
    cache === nothing || return cache

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

    n_lon = length(lon)
    n_lat = length(lat)
    n_class = size(frac, 3)
    lon_min = Float64(lon[1])
    lat_min = Float64(lat[1])
    lon_step = Float64(lon[2]) - lon_min
    lat_step = Float64(lat[2]) - lat_min

    cache = (
        lon = lon, lat = lat, frac = frac,
        lon_min = lon_min, lon_step = lon_step,
        lat_min = lat_min, lat_step = lat_step,
        n_lon = n_lon, n_lat = n_lat, n_class = n_class,
    )
    _LANDUSE_CACHE[] = cache
    return cache
end

"""
    reload_landuse_fractions!()

Invalidate the in-memory cache. The next `landuse_frac_at` call re-reads
the NetCDF from disk. Useful after regenerating the file with
`preprocess_landuse.jl`.
"""
function reload_landuse_fractions!()
    _LANDUSE_CACHE[] = nothing
    return nothing
end

"""
    set_landuse_data_path!(path)

Override the NetCDF path used by `landuse_frac_at`. Invalidates the
cache. Primarily used by tests to point at a synthetic fixture.
"""
function set_landuse_data_path!(path::AbstractString)
    _LANDUSE_NC_PATH[] = String(path)
    _LANDUSE_CACHE[] = nothing
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
function landuse_frac_at(lon_rad, lat_rad, class_index)
    data = _load_landuse()
    lon_deg = rad2deg(lon_rad)
    lat_deg = rad2deg(lat_rad)
    i_lon = clamp(round(Int, (lon_deg - data.lon_min) / data.lon_step) + 1, 1, data.n_lon)
    i_lat = clamp(round(Int, (lat_deg - data.lat_min) / data.lat_step) + 1, 1, data.n_lat)
    return Float64(data.frac[i_lon, i_lat, class_index])
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
# Called from `couple2(::DryDepositionGasFractionalCoupler, ::GEOSFPCoupler)`
# as `season_at(gp.t_ref + t)`. For Southern-Hemisphere or tropical work,
# override `season` manually instead of binding to this function.
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
