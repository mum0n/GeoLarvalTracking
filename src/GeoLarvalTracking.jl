module GeoLarvalTracking

# ==============================================================================
#  GeoLarvalTracking
#
#  Phase 6 of the ParticleTracking split (SPLIT_PLAN.md 1). This package owns the Lagrangian
#  half: larval biology, DVM, molting, mortality, settlement, depth-stratified tessellation,
#  larval storage, and the larval figures. It reads the ocean by LOCATION and by METADATA
#  from the handover `GeoOceanModel.SimulationResult` describes, and imports nothing from
#  GeoOceanModel -- there is no shared TOML and no import edge, by the contract in 1.
#
#  Everything here was moved from `ParticleTracking/src` unchanged, apart from the module name
#  and the paths in `include`. Nothing was renamed.
#
#  The two packages share nothing but GeoData. Measured cross-domain call edges between the
#  biology half and the ocean half of the original source: zero (SPLIT_PLAN.md 7d).
# ==============================================================================

using
    Random,
    Statistics,
    LinearAlgebra,
    TOML,
    JLD2,
    Interpolations,
    Dates,
    Downloads,
    GeoData

# GeoDataSources supplies the data layer for reading a coastline and for a bathymetry
# product, which a seed-position filter needs to reject land. Nothing here fetches over the
# network except through the layer that owns the network.
using GeoDataSources: available_data, DataSource

# The provenance registry, as a submodule. The larval half records what it read for its own
# output, so it carries the same registry the ocean half does -- and this is the one file the
# two packages are permitted to have in common, because it is the data layer's, not either
# half's.
include("data/manifest.jl")
using .Manifest

using Oceananigans.Units
using Oceananigans.Utils: prettytime
using Oceananigans

# ============================================================================
# Source layout, in dependency order: config first because everything falls back on it, the
# flow reader next because the tracking half cannot do anything without a velocity field, and
# then the biology and analysis halves.
# ============================================================================
include("config/configuration.jl")
include("data/store_reader.jl")
include("driver.jl")
include("biology/larval_behavior.jl")
include("analysis/voronoi_tessellation.jl")
include("analysis/larval_storage.jl")
include("analysis/empirical_analysis.jl")
include("output/visualization.jl")
include("output/biology_diagnostics.jl")

end # module
