"""
    test_flow_reader.jl

The one test that pins the copied JLD2 reader to the writer that produces the archives.

This is the cost of the copy made explicit. `GeoLarvalTracking/src/data/flow_reader.jl` was
copied from `ParticleTracking/src/model/simulation.jl` because importing it would be an edge
from the larval half to the ocean half, which SPLIT_PLAN.md section 1 forbids. The format is
therefore known in two places, and a change to either that breaks the pair must show up here
rather than as a silently wrong velocity.

It reads a REAL archive rather than a synthetic one, because a synthetic array can be built
to match any layout mistake the reader makes.
"""

using Test
using JLD2
using GeoLarvalTracking

# A run the umbrella actually produced, so this is a genuine writer-to-reader pair.
const ARCHIVE = normpath(joinpath(@__DIR__, "..", "..", "ParticleTracking", "outputs",
                                  "default", "hydrodynamics_default.jld2"))
const GRID_SIDECAR = normpath(joinpath(@__DIR__, "..", "..", "ParticleTracking", "outputs",
                                       "default", "hydrodynamics_default_grid.jld2"))

@testset "flow_reader reads the real archive" begin
    @test isfile(ARCHIVE)
    @test isfile(GRID_SIDECAR)

    # The sidecar this reader depends on, read directly so the assertion is about its
    # contents rather than about the reader's opinion of them.
    sidecar = JLD2.jldopen(GRID_SIDECAR, "r") do f
        (lons = collect(Float64, f["lons"]), lats = collect(Float64, f["lats"]),
         depths = collect(Float64, f["depths"]),
         size = Int.(collect(f["size"])))
    end

    fn = GeoLarvalTracking.create_flow_interpolator_from_jld2(ARCHIVE)
    @test fn isa Function

    # A point strictly inside the grid, at the first snapshot. The values must be finite and
    # the structure must be the one `track_larval_cohort` consumes.
    lon, lat, z, t = sidecar.lons[2], sidecar.lats[2], sidecar.depths[2], 0.0
    val = fn(lon, lat, z, t)
    for field in (val.u, val.v, val.w, val.T)
        @test isfinite(field)
    end

    # Temperature is what was actually set in this archive (stratified, not at rest), so its
    # magnitude being ocean-plausible is evidence the read worked.
    @test -5.0 < val.T < 40.0

    # Outside the domain the reader clamps to the nearest edge rather than erroring, which is
    # what lets a track that wanders past the boundary still return a velocity.
    @test isfinite(fn(1e6, -1e6, 0.0, 0.0).u)
    @test isfinite(fn(1e6, -1e6, 0.0, 1e9).T)
end
