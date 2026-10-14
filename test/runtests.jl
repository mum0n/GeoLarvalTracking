"""
    runtests.jl -- the larval half's suite

Covers what `GeoLarvalTracking` owns: the store reader that replaces the JLD2-plus-sidecar
reader, and the driver segments moved out of the umbrella script.

The store reader is the coupling with `GeoOceanModel`, and it is tested against a store the
ocean half actually wrote. `GeoLarvalTracking` cannot import `GeoOceanModel` -- SPLIT_PLAN.md
section 1 forbids an edge in that direction -- so the fixture is committed under `test/fixtures/`
by `GeoOceanModel` itself, and this suite reads it. Regenerate it with `GeoOceanModel` rather
than editing it by hand.

A synthetic store is deliberately not used: it can be built to match any layout mistake the
reader makes, so reading the real writer's output is what pins the two halves together.
"""

using Test
using GeoLarvalTracking
using GeoData
using JSON3

const FIXTURE = normpath(joinpath(@__DIR__, "fixtures", "hydrodynamics.zarr",
                                  "hydrodynamics.zarr"))

@testset "GeoLarvalTracking" begin

    @testset "store reader: the ocean half's output, read by this half" begin
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
    end

    @testset "store reader: absence and incomplete metadata are errors" begin
        @test_throws Exception GeoLarvalTracking.open_hydrodynamics(
            normpath(joinpath(@__DIR__, "fixtures", "absent.zarr")))
        # A field the store does not carry, rather than zeros passed off as a velocity.
        @test_throws Exception GeoLarvalTracking.create_flow_interpolator(
            FIXTURE; fields = ("u", "nonsense"))
        # An interval beyond the run's length, rather than the last interval.
        local s = GeoLarvalTracking.open_hydrodynamics(FIXTURE)
        @test_throws Exception GeoLarvalTracking.variable(s, "T", 10_000)
    end

    @testset "driver: the segments are reachable from the package" begin
        # These lived in the umbrella's script, where `Pkg.test` could not reach them. Being
        # defined in `src/` is what makes them testable at all.
        for name in ("run_segment_tracking", "run_segment_metrics", "run_segment_visualize",
                     "render_all_figures", "load_run_trajectories", "run_cli_list_runs",
                     "run_cli_compare_scenarios", "run_cli_model_average",
                     "create_flow_interpolator", "open_hydrodynamics")
            @test isdefined(GeoLarvalTracking, Symbol(name))
        end
    end
end
