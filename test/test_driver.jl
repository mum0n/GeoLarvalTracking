using Test
using GeoLarvalTracking

@testset "Driver Segments & Surface API" begin
    # These lived in the umbrella's script, where `Pkg.test` could not reach them. Being
    # defined in `src/` is what makes them testable at all.
    for name in ("run_segment_tracking", "run_segment_metrics", "run_segment_visualize",
                 "render_all_figures", "load_run_trajectories", "run_cli_list_runs",
                 "run_cli_compare_runs", "run_cli_model_average",
                 "create_flow_interpolator", "open_hydrodynamics")
        @test isdefined(GeoLarvalTracking, Symbol(name))
    end
end
