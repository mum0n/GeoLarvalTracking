using Test
using GeoLarvalTracking
using GeoData

const FIXTURE = normpath(joinpath(@__DIR__, "fixtures", "hydrodynamics.zarr",
                                  "hydrodynamics.zarr"))

@testset "Store Reader" begin
    @test isdir(FIXTURE)
    local s = GeoLarvalTracking.open_hydrodynamics(FIXTURE)
    @test sort(collect(s.variables)) == ["T", "u", "v", "w"]
    @test length(s.axes.lon) == 5
    @test length(s.axes.lat) == 4
    @test length(s.axes.depth) == 2

    # The time axis is the store's own coordinate, read back through GeoData rather than
    # assumed -- and rather than from a `time.json` sidecar, which the writer no longer
    # keeps now that the axis lives inside the dataset.
    local back = GeoData.geoload(FIXTURE; backend = :zarr)
    @test s.times == collect(Float64, back.coords[:time].data)
    @test length(s.times) == 3

    # A field at an interval is the value written.
    local T0 = GeoLarvalTracking.variable(s, "T", 1)
    @test size(T0) == (5, 4, 2)
    @test T0[1, 1, 1] == 12.0

    # The interpolator the tracker integrates against.
    local fn = GeoLarvalTracking.create_flow_interpolator(FIXTURE)
    local v = fn(0.0, 0.0, -2.0, 0.0)
    @test v.T == 12.0
    @test isfinite(v.u) && isfinite(v.v) && isfinite(v.w)

    # Outside the domain it clamps rather than erroring, so a track that wanders past the
    # boundary still returns a real velocity.
    @test isfinite(fn(1e6, -1e6, 0.0, 0.0).T)

    # Absence and incomplete metadata are errors
    @test_throws Exception GeoLarvalTracking.open_hydrodynamics(
        normpath(joinpath(@__DIR__, "fixtures", "absent.zarr")))
    @test_throws Exception GeoLarvalTracking.create_flow_interpolator(
        FIXTURE; fields = ("u", "nonsense"))
    local s2 = GeoLarvalTracking.open_hydrodynamics(FIXTURE)
    @test_throws Exception GeoLarvalTracking.variable(s2, "T", 10_000)
end
