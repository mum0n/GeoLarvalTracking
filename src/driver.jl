using Dates
using Random
using Statistics

# `Build` is a local helper of the umbrella driver, imported unqualified by one of the moved
# functions. It is not this package's, so it is stubbed rather than left to resolve at call time.
Build(args...; kwargs...) = error("`Build` is a local helper of the umbrella driver, not this " *
                                  "package's: the extraction carried its docstring without its " *
                                  "body. Whoever needs it should port it deliberately.")
    """

resolve_hydro_model_path(opts::HydrodynamicOptions, default_filename::String) -> Tuple{String, String}

Resolve the target hydrodynamic JLD2 output/input path and filename. If `opts.hydro_model_file`
is specified, it is used directly (or resolved relative to `opts.output_dir` if a bare filename
is provided). If empty, returns `(joinpath(opts.output_dir, default_filename), default_filename)`.

# Inputs
- `opts::HydrodynamicOptions`: Configuration parameters.
- `default_filename::String`: Default fallback filename when none is specified.

# Outputs
- `Tuple{String, String}`: `(full_path, file_basename)`
"""
function resolve_hydro_model_path(opts::HydrodynamicOptions, default_filename::String)
    if isempty(opts.hydro_model_file)
        full_path = joinpath(opts.output_dir, default_filename)
        return (full_path, default_filename)
    else
        raw = opts.hydro_model_file
        full_path = isabspath(raw) || dirname(raw) != "" ? raw : joinpath(opts.output_dir, raw)
        return (full_path, basename(raw))
    end
end


"""
    run_segment_tracking(; opts::HydrodynamicOptions) -> NamedTuple

Simulate Lagrangian particle transport for a snow crab larval cohort across nursery
spawning areas with DVM swimming, M2 tidal currents, degree-day molting, and
benthic settlement filtering.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.

# Outputs
- `NamedTuple`: Trajectories record and summary metrics.
"""
function run_segment_tracking(; opts::HydrodynamicOptions = HydrodynamicOptions())
    println("\n=================================================================")
    println(" [Segment 6/8] Lagrangian Particle Tracking & Larval Behavior")
    println("=================================================================")
    rng = MersenneTwister(opts.seed)

    # Spawning area: derive from loaded CFA boundaries + user-defined spatial buffer
    cfa_polys = load_cfa_polygons(opts.input_dir)
    spawn_lon, spawn_lat = if !isempty(cfa_polys)
        env = get_strata_buffered_envelope(cfa_polys, buffer_km = opts.buffer_km)
        (env.lon_range, env.lat_range)
    else
        expand_domain_with_buffer(opts.domain_lon, opts.domain_lat, buffer_km = opts.buffer_km)
    end

    # Active seabed bathymetry dataset. The run's resolved sea floor, as a `(lon, lat) -> z`
    # function, rather than a path: which file holds it is `GeoDataSources`'s decision.
    bathy_src = bathy_field_fn

    println("Initializing $(opts.n_particles) Zoea I larvae in marine water (mode = $(opts.release_depth_mode), min depth >= $(opts.min_seabed_depth) m)...")
    larvae = initialize_larval_particles(
        opts.n_particles,
        lon_range = spawn_lon,
        lat_range = spawn_lat,
        release_depth_mode = opts.release_depth_mode,
        bottom_offset = opts.bottom_release_offset,
        ascent_target_depth = opts.ascent_target_depth,
        min_seabed_depth = opts.min_seabed_depth,
        buffer_km = 0.0, # Buffer is already incorporated into spawn_lon/spawn_lat
        bathymetry = bathy_src,
        stage = :zoea1,
        rng = rng
    )

    # 4D hydrodynamic and temperature fields: use create_flow_interpolator_from_jld2
    # if simulation JLD2 archive exists; otherwise use analytical Scotian Shelf background.
    default_jld2 = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
    jld2_target_path, _ = resolve_hydro_model_path(opts, default_jld2)

    if opts.track_only && !isfile(jld2_target_path)
        error("Cannot run in --track-only mode: hydrodynamic model file does not exist: $(jld2_target_path)\n" *
              "Please run with --hydro-only first or specify an existing file with --hydro-model=<path>.")
    end

    sim_jld2_candidates = unique([
        jld2_target_path,
        joinpath(opts.output_dir, default_jld2),
        joinpath(opts.output_dir, "simulation_flow.jld2"),
        joinpath(opts.output_dir, "simulation_$(opts.scenario)_$(opts.projection_year).jld2"),
        joinpath("outputs", "simulation_flow.jld2")
    ])
    sim_jld2_path = findfirst(isfile, sim_jld2_candidates)
    # Span, in seconds, of the hydrodynamic record the cohort will be advected through. The cohort
    # is tracked for this whole span and stops early only once every larva is dead or settled, so
    # the drift is never shorter than the model it is driven by.
    hydro_span = [0.0, 0.0]

    flow_interpolator = if !isnothing(sim_jld2_path)
        actual_path = sim_jld2_candidates[sim_jld2_path]
        println("Loading 4D simulated hydrodynamic fields from $(actual_path)...")
        try
            # A production run stores thousands of snapshots; materialising all of them as Float64
            # for u, v, w and T is hundreds of gigabytes and exhausts memory, which previously made
            # this silently fall through to the analytical jet below -- producing recruitment numbers
            # from a synthetic flow rather than the simulated one. So: keep the full record's time
            # span, but bound how many snapshots are held resident via `max_flow_snapshots`
            # (evenly subsampled, endpoints preserved), and report the retained span back.
            interp = create_flow_interpolator(actual_path;
                                                        t_window = nothing,
                                                        max_snapshots = opts.max_flow_snapshots,
                                                        time_range_out = hydro_span)
            println("  retained span: $(round(hydro_span[1], digits=1)) .. " *
                    "$(round(hydro_span[2], digits=1)) s " *
                    "($(round((hydro_span[2] - hydro_span[1]) / 86400, digits=2)) days)")
            interp
        catch err
            # Falling back to an analytical current here would silently replace the simulated
            # hydrodynamics with a synthetic jet, turning a memory limit into a wrong scientific
            # result. Make the failure explicit and actionable instead.
            if opts.allow_analytical_fallback
                @warn "Failed to build the flow interpolator from $(actual_path): " *
                      "$(sprint(showerror, err)). Falling back to the analytical jet because " *
                      "allow_analytical_fallback is set -- results will NOT reflect the simulation."
                nothing
            else
                error(
                    "Failed to build the flow interpolator from $(actual_path):\n  " *
                    sprint(showerror, err) * "\n\n" *
                    "Tracking needs the simulated fields; substituting an analytical current would " *
                    "silently change the science. To proceed deliberately with a synthetic flow, " *
                    "re-run with --allow-analytical-fallback. To reduce memory, lower " *
                    "[hydrodynamics].max_flow_snapshots (currently $(opts.max_flow_snapshots))."
                )
            end
        end
    else
        nothing
    end

    flow_field_fn = if !isnothing(flow_interpolator)
        (lon, lat, z, t) -> begin
            f = flow_interpolator(lon, lat, z, t)
            (f.u, f.v, f.w)
        end
    else
        # Scotian Shelf alongshore current: flows southwestward along the shelf edge.
        # u ≈ -0.08 to -0.12 m/s (westward), v ≈ -0.03 to -0.06 m/s (southward/offshore).
        # Tidal oscillation superimposed with ~6 h period at quarter-amplitude.
        # Reference: Loder, J. W. & Petrie, B. (1991), CJFAS.
        (lon, lat, z, t) -> begin
            tidal_phase = 2.0 * π * t / 44712.0  # M2 period ≈ 12.42 h
            depth_decay  = exp(z / 120.0)          # velocity decays with depth
            lon_norm = clamp((lon - opts.domain_lon[1]) /
                             (opts.domain_lon[2] - opts.domain_lon[1]), 0.0, 1.0)
            jet_factor = 0.5 + 0.8 * exp(-((lon_norm - 0.3)^2) / 0.04)
            u_mean = (-0.10 * jet_factor + 0.02 * sin(tidal_phase)) * depth_decay
            v_mean = (-0.04 * jet_factor + 0.01 * cos(tidal_phase)) * depth_decay
            w_mean = 0.00020 * sin(2.0 * π * lon_norm) * depth_decay
            (u_mean, v_mean, w_mean)
        end
    end

    temp_field_fn = if !isnothing(flow_interpolator)
        (lon, lat, z, t) -> begin
            f = flow_interpolator(lon, lat, z, t)
            f.T
        end
    else
        # 3D Scotian Shelf thermal structure: surface mixed layer + CIL + slope water.
        # Consistent with set_initial_stratification! in hydrodynamic_model.jl.
        (lon, lat, z, t) -> begin
            lon_min, lon_max = opts.domain_lon
            lat_min, lat_max = opts.domain_lat
            x_norm = clamp((lon - lon_min) / (lon_max - lon_min), 0.0, 1.0)
            y_norm = clamp((lat - lat_min) / (lat_max - lat_min), 0.0, 1.0)
            T_surf   = 14.0 + 8.0 * x_norm - 5.0 * y_norm   # surface mixed layer
            T_cil_min = 1.5
            T_slope   = 8.5
            if z > -20.0
                frac = (z + 20.0) / 20.0
                clamp(T_cil_min + frac * (T_surf - T_cil_min), -1.5, 22.0)
            elseif z > -80.0
                centre_frac = (z + 50.0) / 30.0
                clamp(T_cil_min + 2.5 * centre_frac^2, -1.5, T_surf)
            else
                depth_factor = (abs(z) - 80.0) / 100.0
                clamp(T_cil_min + depth_factor * (T_slope - T_cil_min), T_cil_min, T_slope)
            end
        end
    end

    # Seabed bathymetric elevation. Resolved from the run's source and range through
    # `GeoDataSources`, so the seeded larvae sit on the same sea floor the model was built on.
    # The interpolator is built from those arrays rather than from a file, so nothing here
    # opens a NetCDF, and a run that never fetched one gets a real product rather than an
    # analytic stand-in.
    run_bathy = _run_bathymetry(; opts = opts)
    bathy_field_fn = get_bathymetry_interpolator(
        run_bathy.elevation, Float64.(collect(run_bathy.lon)), Float64.(collect(run_bathy.lat)))

    coastline_path = joinpath(opts.input_dir, "coastline.dat")
    coast_polys = isfile(coastline_path) ? load_coastline_polygons(coastline_path) : nothing

    # Track over the full span of the hydrodynamic record. The cohort stops early by itself once
    # every larva is dead or settled, so this is an upper bound the run only reaches if larvae
    # persist. When there is no hydro file (analytical fallback or synthetic quick run) the
    # configured pelagic larval duration is the only sensible horizon.
    track_horizon = if !isnothing(flow_interpolator) && hydro_span[2] > hydro_span[1]
        hydro_span[2] - hydro_span[1]
    else
        opts.track_duration
    end
    println("Tracking cohort over up to $(track_horizon / 86400.0) days " *
            "(dt=$(opts.track_dt)s, stops early when all larvae are dead or settled)...")
    trajectories = track_larval_cohort(
        larvae,
        velocity_fn = flow_field_fn,
        temperature_fn = temp_field_fn,
        bathymetry_fn = bathy_field_fn,
        total_duration = track_horizon,
        dt = opts.track_dt,
        max_current_speed = opts.max_current_speed,
        κ_h = opts.diffusivity_h,
        κ_v = opts.diffusivity_v,
        is_lat_lon = true,
        enable_tides = opts.enable_tides,
        tidal_u_amp = opts.tidal_u_amp,
        tidal_v_amp = opts.tidal_v_amp,
        # The constituents and their amplitude ratios come from the TOML via `opts`. Hardcoding
        # them here would run the tides on a set the operator did not choose, which is the
        # same class of silent substitution the data path was changed to avoid.
        tidal_constituents = opts.tidal_constituents,
        tidal_u_amplitudes = Dict(:M2 => opts.tidal_u_amp,
                                  :S2 => 0.44 * opts.tidal_u_amp,
                                  :N2 => 0.19 * opts.tidal_u_amp,
                                  :K1 => 0.18 * opts.tidal_u_amp,
                                  :O1 => 0.15 * opts.tidal_u_amp),
        tidal_v_amplitudes = Dict(:M2 => opts.tidal_v_amp,
                                  :S2 => 0.42 * opts.tidal_v_amp,
                                  :N2 => 0.19 * opts.tidal_v_amp,
                                  :K1 => 0.18 * opts.tidal_v_amp,
                                  :O1 => 0.15 * opts.tidal_v_amp),
        enable_molting = opts.enable_molting,
        # Stochasticity controls. These must be passed explicitly: `track_larval_cohort` defaults
        # all three coefficients of variation to 0.0, i.e. a fully deterministic cohort. Loading
        # them into `opts` and writing them to `resolved_config.toml` is not enough on its own --
        # if they are not forwarded here the provenance file records a stochastic model that the
        # run never actually used.
        cv_molt = opts.cv_molt,
        cv_mortality = opts.cv_mortality,
        cv_settlement = opts.cv_settlement,
        settlement_stochastic = opts.settlement_stochastic,
        enable_bbl = true,
        enable_sinking = true,
        enable_initial_ascent = opts.enable_initial_ascent,
        ascent_speed = opts.ascent_speed,
        ascent_target_depth = opts.ascent_target_depth,
        coastline = coast_polys,
        rng = rng
    )

    # Save trajectories to JLD2 checkpoint for modular reloading
    mkpath(opts.output_dir)
    tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
    track_checkpoint = joinpath(opts.output_dir, "larval_trajectories$(tag).jld2")
    jldsave(
        track_checkpoint;
        lons = trajectories.lons,
        lats = trajectories.lats,
        depths = trajectories.depths,
        temperatures = hasproperty(trajectories, :temperatures) ? trajectories.temperatures : nothing,
        degree_days = trajectories.degree_days,
        degree_days_timeseries = hasproperty(trajectories, :degree_days_timeseries) ? trajectories.degree_days_timeseries : nothing,
        survival_probability = hasproperty(trajectories, :survival_probability) ? trajectories.survival_probability : nothing,
        stage_survival = hasproperty(trajectories, :stage_survival) ? trajectories.stage_survival : nothing,
        stages = trajectories.stages,
        alive = trajectories.alive,
        settlement_status = trajectories.settlement_status,
        settlement_age = trajectories.settlement_age,
        ascent_duration = hasproperty(trajectories, :ascent_duration) ? trajectories.ascent_duration : nothing,
        # The cohort-level stage diagnostics are saved too, so `--viz` can redraw the
        # progression figures from a checkpoint without re-running Segment 6.
        cohort_molt_fraction = hasproperty(trajectories, :cohort_molt_fraction) ? trajectories.cohort_molt_fraction : nothing,
        molt_schedule = hasproperty(trajectories, :molt_schedule) ? trajectories.molt_schedule : nothing,
        times = trajectories.times,
        ids = trajectories.ids
    )

    # Maintain default checkpoint as fallback when custom run_id is supplied
    if !isempty(opts.run_id)
        default_cp = joinpath(opts.output_dir, "larval_trajectories.jld2")
        try
            cp(track_checkpoint, default_cp, force = true)
        catch
        end
    end

    # Save trajectories directly into Zarr
    if opts.enable_zarr
        try
            println("Archiving trajectories to Zarr -> $(opts.zarr_path)...")
            db = open_storage(opts.zarr_path)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                save_simulation_run!(
                    db, target_run_id, opts;
                    trajectories = trajectories,
                    config = options_to_configuration(opts),
                    notes = "Hydrodynamic tracking run ($(opts.scenario), $(opts.projection_year))" *
                            (!isempty(opts.run_id) ? " [$(opts.run_id)]" : "")
                )
                println("Trajectories for '$(target_run_id)' successfully archived in Zarr.")
            finally
                close_storage(db)
            end
        catch err
            @warn "Failed to archive trajectories to Zarr: $(err)"
        end
    end

    n_settled = count(==( :settled_successful), trajectories.settlement_status)
    n_alive   = count(identity, trajectories.alive)
    println("Lagrangian tracking complete:")
    println("  Total particles:     $(opts.n_particles)")
    println("  Surviving particles: $(n_alive) / $(opts.n_particles)")
    println("  Settled on nursery:  $(n_settled) / $(opts.n_particles) " *
            "($(round(100.0 * n_settled / opts.n_particles, digits=1))%)")
    println("  Checkpoint saved:    $(track_checkpoint)")

    return (trajectories = trajectories, checkpoint_path = track_checkpoint)
end


"""
    run_segment_metrics(;
        opts::HydrodynamicOptions,
        trajectories=nothing
    ) -> NamedTuple

Compute gridded larval retention and recruitment metrics, thermal exposure indices,
empirical Lagrangian drift velocities and effective diffusivities, regional
connectivity matrices across Crab Fishing Areas (CFAs), and export multi-layer NetCDF and JLD2.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.
- `trajectories`: Optional trajectory NamedTuple.

# Outputs
- `NamedTuple`: Metric outputs, connectivity matrix, and export filepaths.
"""
function run_segment_metrics(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    trajectories = nothing
)
    target_trajs = isnothing(trajectories) ? load_run_trajectories(opts) : trajectories

    println("\n=================================================================")
    println(" [Segment 7/8] Empirical Movement, Recruitment & Connectivity")
    println("=================================================================")

    lon_b = range(opts.domain_lon[1], opts.domain_lon[2], length = opts.grid_size[1])
    lat_b = range(opts.domain_lat[1], opts.domain_lat[2], length = opts.grid_size[2])

    # 1. Empirical Advection & Turbulent Diffusivity
    println("Estimating empirical advection and turbulent diffusivity fields...")
    emp_mov = estimate_empirical_movement(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 2. Gridded Recruitment & Settlement Metrics
    println("Computing gridded recruitment and benthic nursery settlement metrics...")
    rec_metrics = compute_gridded_recruitment_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 3. Gridded Thermal Exposure & Degree-Days
    println("Computing gridded thermal exposure metrics...")
    therm_metrics = compute_gridded_thermal_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 4. Demographic Connectivity Matrix across Crab Fishing Areas (CFAs)
    println("Computing macro-regional demographic connectivity matrix...")
    cfa_polys = load_cfa_polygons(opts.input_dir)
    cfa_definitions = if !isempty(cfa_polys)
        println("Loaded $(length(cfa_polys)) CFA boundary polygons from $(opts.input_dir)/: $(join([p.name for p in cfa_polys], ", "))")
        vcat(cfa_polys, [(name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))])
    else
        [
            (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
            (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
            (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
            (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
        ]
    end
    conn = compute_empirical_connectivity(target_trajs, strata_definitions = cfa_definitions)
    println("Transition Probability Matrix (P_ij):")
    for i in 1:length(conn.strata_names)
        row_str = join([string(round(conn.matrix[i, j], digits = 3)) for j in 1:length(conn.strata_names)], ", ")
        println("  $(rpad(conn.strata_names[i], 28)) -> [$(row_str)]")
    end

    # 5. Depth-stratified Voronoi Tessellation & Multi-Resolution Connectivity
    voronoi_res = if opts.enable_voronoi
        println("\nComputing depth-stratified Voronoi tessellation ($(opts.voronoi_n_units) units)...")
        v_tess = generate_depth_stratified_voronoi_units(
            :numerical_earth;
            lon_range = opts.domain_lon,
            lat_range = opts.domain_lat,
            n_units = opts.voronoi_n_units,
            prob_core = opts.voronoi_prob_core,
            prob_shallow = opts.voronoi_prob_shallow,
            prob_deep = opts.voronoi_prob_deep,
            min_res_core_km = opts.voronoi_min_res_core_km,
            min_res_shallow_km = opts.voronoi_min_res_shallow_km,
            min_res_deep_km = opts.voronoi_min_res_deep_km,
            # The TOML states band edges as positive depths below the surface; the tessellator
            # bins on signed z, so the deeper edge gets the sign.
            core_depth_range = (-opts.voronoi_core_depth_max, -opts.voronoi_core_depth_min),
            shallow_depth_range = (-opts.voronoi_shallow_depth_max, -opts.voronoi_shallow_depth_min),
            deep_depth_range = (-opts.voronoi_deep_depth_max, -opts.voronoi_deep_depth_min),
            slope_weighting = opts.voronoi_slope_weighting,
            slope_factor = opts.voronoi_slope_factor,
            seed = opts.seed
        )
        println("Generated $(length(v_tess.units)) Voronoi units across depth strata.")

        # Export Voronoi polygons to standard GeoJSON format for GIS integration
        geojson_export_path = joinpath(opts.output_dir, "voronoi_units.geojson")
        try
            export_voronoi_geojson(v_tess, geojson_export_path)
            println("Exported Voronoi polygons to GeoJSON -> $(geojson_export_path)")
        catch geo_err
            @warn "Failed exporting Voronoi GeoJSON: $(geo_err)"
        end

        v_metrics = compute_tesselated_connectivity_matrix(target_trajs, v_tess)
        println("Voronoi Macro-Strata Transition Matrix (3x3):")
        for (idx, st) in enumerate(v_metrics.strata_names)
            row_str = join([string(round(v_metrics.macro_matrix[idx, j], digits = 3)) for j in 1:length(v_metrics.strata_names)], ", ")
            println("  $(rpad(string(st), 12)) -> [$(row_str)]")
        end
        (tessellation = v_tess, metrics = v_metrics)
    else
        nothing
    end

    # 6. Export comprehensive multi-variable NetCDF and JLD2 archives
    nc_export_path = joinpath(opts.output_dir, "larval_dispersal_analysis.nc")
    jld_export_path = joinpath(opts.output_dir, "larval_dispersal_analysis.jld2")
    active_config = options_to_configuration(opts)

    println("Exporting multi-layer NetCDF archive -> $(nc_export_path)...")
    export_larval_dispersal_netcdf(
        nc_export_path,
        trajectories = target_trajs,
        lon_bins = lon_b,
        lat_bins = lat_b,
        strata_definitions = cfa_definitions,
        config = active_config
    )

    println("Exporting comprehensive JLD2 archive -> $(jld_export_path)...")
    export_larval_dispersal_jld2(
        jld_export_path,
        trajectories = target_trajs,
        lon_bins = lon_b,
        lat_bins = lat_b,
        strata_definitions = cfa_definitions,
        config = active_config
    )

# 7. Archive simulation run, trajectories, metrics & connectivity in Zarr
    if opts.enable_zarr
        try
            println("Archiving simulation run and metrics to Zarr -> $(opts.zarr_path)...")
            db = open_storage(opts.zarr_path)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                save_simulation_run!(
                    db,
                    target_run_id,
                    opts;
                    trajectories = target_trajs,
                    metrics = (
                        mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                        mean_degree_days = therm_metrics.mean_degree_days
                    ),
                    connectivity = conn,
                    gridded_dispersal = (
                        lon_centers = emp_mov.lon_centers,
                        lat_centers = emp_mov.lat_centers,
                        u_mean = emp_mov.u_mean,
                        v_mean = emp_mov.v_mean,
                        diffusivity = emp_mov.diffusivity,
                        density = rec_metrics.settlement_density,
                        mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                        mean_degree_days = therm_metrics.mean_degree_days,
                        sample_count = emp_mov.sample_count
                    ),
                    config = active_config,
                    notes = "Hydrodynamic workflow run ($(opts.scenario), $(opts.projection_year))" *
                            (!isempty(opts.run_id) ? " [$(opts.run_id)]" : "")
                )
                if opts.enable_voronoi && !isnothing(voronoi_res)
                    try
                        # Zarr doesn't support raw SQL - use the storage API instead
                        if !haskey(db["metrics"], target_run_id)
                            Zarr.create_group(db["metrics"], target_run_id)
                        end
                        mgroup = db["metrics"][target_run_id]
                        for u in voronoi_res.units
                            row_id = u.id
                            Zarr.write(mgroup["unit_$(row_id)"], Dict(
                                "run_id" => target_run_id,
                                "unit_id" => u.id,
                                "lon" => u.lon,
                                "lat" => u.lat,
                                "depth" => u.depth,
                                "stratum" => string(u.stratum),
                                "area_km2" => u.area_km2,
                                "settlement_count" => u.settlement_count,
                            ))
                        end
                        println("Archived $(length(voronoi_res.units)) Voronoi units to Zarr.")
                    catch v_err
                        @warn "Failed to archive Voronoi units to Zarr: $(v_err)"
                    end
                end
                println("Zarr run '$(target_run_id)' successfully archived.")
            finally
                close_storage(db)
            end
        catch err
            @warn "Failed to archive simulation run and metrics to Zarr: $(err)"
        end
    end

    return (
        empirical_movement = emp_mov,
        recruitment_metrics = rec_metrics,
        thermal_metrics = therm_metrics,
        connectivity = conn,
        voronoi = voronoi_res,
        netcdf_path = nc_export_path,
        jld2_path = jld_export_path,
        storage_path = opts.enable_zarr ? opts.zarr_path : nothing
    )
end


"""
    render_all_figures(; opts) -> NamedTuple

Regenerate every figure and visual for an existing run, from stored results only.

This is the implementation behind `--figures`: particle figures, the larval-biology diagnostics,
the Eulerian hydrodynamic set, the interactive map, and the animation. It reads the trajectory
checkpoint and the hydrodynamic archive rather than re-simulating, so it is the fast path to
iterate on a plotting function.

# What it draws
- **Larval** (7): trajectories, DVM profiles, settlement density, empirical movement field,
  connectivity matrix, thermal exposure, recruitment summary, plus the three biology diagnostics
  (molt progression, degree-day growth, survival curves).
- **Eulerian** (9): advection, tracers, stratification, diffusion, section, station time series,
  time–depth diagram, near-bed temperature map, T–S diagram.
- **Interactive**: the decimated HTML dashboard.
- **Animation**: the field/dashboard animation, when `--animate-hydro` is set.

Missing inputs are reported rather than fatal, so a run with, say, no hydrodynamic archive still
produces the larval figures.
"""
function render_all_figures(; opts::HydrodynamicOptions = HydrodynamicOptions())
    println("\n=================================================================")
    println(" [Figures] Regenerating all figures and visuals")
    println("=================================================================")
    viz = run_segment_visualize(opts = opts)
    if opts.animate_hydro
        anim_path = if !isempty(opts.hydro_model_file) && isfile(opts.hydro_model_file)
            render_hydrodynamic_animation_if_requested(opts, opts.hydro_model_file)
        else
            println("  --animate-hydro set but no hydro model file available; skipping animation.")
            nothing
        end
        return merge(viz, (animation = anim_path,))
    end
    return viz
end


"""
    _ck_or(saved::Dict, key::AbstractString, default)

Read `key` from a loaded JLD2 checkpoint dictionary, falling back to `default` when the key is
absent **or** stored as `nothing`.

The tracking checkpoint writes optional fields as `nothing` when the trajectory lacks them, and
older checkpoints predate several fields entirely, so a plain `get(saved, key, default)` is not
enough — a stored `nothing` would flow into the figures as a missing value.
"""
_ck_or(saved, key::AbstractString, default) =
    (haskey(saved, key) && !isnothing(saved[key])) ? saved[key] : default


"""
    _trajectories_from_checkpoint(saved::Dict) -> NamedTuple

Build the trajectory NamedTuple from a loaded `larval_trajectories*.jld2` dictionary.

Fields that older checkpoints lack are filled with the same placeholders used before
(`temperatures` 4 °C, `degree_days_timeseries` 40 °C·d, `survival_probability` 0.95). The
stage diagnostics `cohort_molt_fraction` and `molt_schedule` are set to `nothing` when absent
and are then rebuilt by [`_complete_stage_diagnostics`](@ref); a zero placeholder would plot
as a cohort that never moults.
"""
function _trajectories_from_checkpoint(saved)
    n_p, n_t = size(saved["lons"])
    t_end = saved["times"][end]
    return (
        lons = saved["lons"],
        lats = saved["lats"],
        depths = saved["depths"],
        temperatures = _ck_or(saved, "temperatures", fill(4.0, n_p, n_t)),
        degree_days = saved["degree_days"],
        degree_days_timeseries = _ck_or(saved, "degree_days_timeseries",
                                        fill(40.0, n_p, n_t)),
        survival_probability = _ck_or(saved, "survival_probability", fill(0.95, n_p, n_t)),
        stages = saved["stages"],
        alive = saved["alive"],
        settlement_status = saved["settlement_status"],
        settlement_age = _ck_or(saved, "settlement_age", fill(t_end, n_p)),
        stage_survival = _ck_or(saved, "stage_survival", fill(1.0, n_p)),
        ascent_duration = _ck_or(saved, "ascent_duration", fill(0.0, n_p)),
        cohort_molt_fraction = _ck_or(saved, "cohort_molt_fraction", nothing),
        molt_schedule = _ck_or(saved, "molt_schedule", nothing),
        times = saved["times"],
        ids = saved["ids"]
    )
end


"""
    _complete_stage_diagnostics(trajs::NamedTuple, opts) -> NamedTuple

Add `molt_schedule` and `cohort_molt_fraction` to a trajectory set that lacks them, such as one
read back from Zarr/GeoData storage, which stores per-particle rows only.

Both are deterministic functions of the run configuration and the stored degree-days, so they
are recomputed rather than approximated:

- `molt_schedule = molt_schedule(cv_molt = opts.cv_molt)`. Segment 6 calls
  `track_larval_cohort` without degree-day thresholds, so the run used the shared defaults
  (65, 130, 200 °C·d); if Segment 6 starts passing thresholds, pass them here too.
- `cohort_molt_fraction[:, j] = cohort_molt_fraction(schedule, D̄_j)` with
  ``\\bar D_j = N_p^{-1} \\sum_p D_{p,j}``, the cohort-mean degree-day at step `j`. This is the
  formula `track_larval_cohort` applies, and `degree_days_timeseries[:, j]` holds the same
  per-larva values it averages (dead larvae carry their last value).

Fields already present are left unchanged. When molting is disabled the run recorded no
schedule, so nothing is added and the figure reports it as not recorded.
"""
function _complete_stage_diagnostics(trajs::NamedTuple, opts)
    has(k) = hasproperty(trajs, k) && !isnothing(getproperty(trajs, k))
    opts.enable_molting || return trajs

    sched = has(:molt_schedule) ? trajs.molt_schedule : molt_schedule(cv_molt = opts.cv_molt)
    extra = has(:molt_schedule) ? NamedTuple() : (molt_schedule = sched,)

    if !has(:cohort_molt_fraction)
        dds = trajs.degree_days_timeseries
        n_p, n_t = size(dds)
        cm = Matrix{Float64}(undef, 3, n_t)
        for j in 1:n_t
            cm[:, j] .= cohort_molt_fraction(sched, sum(view(dds, :, j)) / n_p)
        end
        extra = merge(extra, (cohort_molt_fraction = cm,))
    end
    return isempty(extra) ? trajs : merge(trajs, extra)
end


"""
    load_run_trajectories(opts::HydrodynamicOptions) -> NamedTuple

Load the trajectory set of the configured run for Segments 7 and 8.

Sources, in order:
1. Analytical storage (`opts.zarr_path`, run id `opts.run_id` or `run_<scenario>_<year>`). If the JLD2
   checkpoint describes the same cohort (equal `ids` and `times`), the checkpoint is used
   instead, because it also holds `stage_survival`, `ascent_duration` and the true
   `settlement_age`, none of which analytical storage stores.
2. The JLD2 checkpoint `larval_trajectories[_<run_id>].jld2`.
3. Neither present: Segment 6 is run.

The result is passed through [`_complete_stage_diagnostics`](@ref), so every source yields the
same set of stage fields.
"""
function load_run_trajectories(opts::HydrodynamicOptions)
    target_run_id = !isempty(opts.run_id) ? opts.run_id :
        "run_$(opts.scenario)_$(opts.projection_year)"
    tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
    tagged_cp = joinpath(opts.output_dir, "larval_trajectories$(tag).jld2")
    track_cp = isfile(tagged_cp) ? tagged_cp :
               joinpath(opts.output_dir, "larval_trajectories.jld2")

    from_db = nothing
    if opts.enable_zarr && isdir(opts.zarr_path)
        db = open_storage(opts.zarr_path; read_only = true)
        try
            println("Loading particle trajectories from Zarr for run '$(target_run_id)'...")
            from_db = load_trajectories_namedtuple(db, target_run_id)
        catch err
            println("  Zarr load failed ($(sprint(showerror, err))); trying the checkpoint.")
        finally
            close_storage(db)
        end
    end

    from_cp = isfile(track_cp) ? _trajectories_from_checkpoint(JLD2.jldopen(track_cp, "r")) : nothing

    trajs = if !isnothing(from_db) && !isnothing(from_cp) &&
               Int.(from_cp.ids) == Int.(from_db.ids) &&
               Float64.(from_cp.times) == Float64.(from_db.times)
        println("  using the matching checkpoint $(track_cp) (holds the full stage record).")
        from_cp
    elseif !isnothing(from_db)
        from_db
    elseif !isnothing(from_cp)
        println("Loading particle trajectories from checkpoint: $(track_cp)...")
        from_cp
    else
        println("Trajectories checkpoint not found. Running Segment 6...")
        run_segment_tracking(opts = opts).trajectories
    end
    return _complete_stage_diagnostics(trajs, opts)
end


"""
    run_segment_visualize(;
        opts::HydrodynamicOptions,
        trajectories=nothing
    ) -> NamedTuple

Generate and export publication-ready CairoMakie figures and standalone
interactive Leaflet HTML dashboards with multi-scenario comparison and
rich spatial layers.

# Outputs
- `NamedTuple`: Paths to generated PNG figures and HTML dashboard.
"""
function run_segment_visualize(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    trajectories = nothing
)
    target_trajs = isnothing(trajectories) ? load_run_trajectories(opts) : trajectories

    println("\n=================================================================")
    println(" [Segment 8/8] Scientific Visualizations & Spatial Figures")
    println("=================================================================")
    mkpath(opts.output_dir)

 

    # Load bathymetry data for background contours. The run's own artefact is used, resolved
    # through `_run_bathymetry`, so the contours and the grid agree about the sea floor rather
    # than one being drawn from whatever file happens to be present.
    bathy_data = _run_bathymetry(; opts = opts)

    lon_b = range(opts.domain_lon[1], opts.domain_lon[2], length = opts.grid_size[1])
    lat_b = range(opts.domain_lat[1], opts.domain_lat[2], length = opts.grid_size[2])

    # 1. 2D Particle Trajectories Map
    fig1_path = joinpath(opts.output_dir, "larval_trajectories.png")
    println("Rendering 2D spatial trajectories map -> $(fig1_path)...")
    cfa_polys = load_cfa_polygons(opts.input_dir)
    plot_particle_trajectories(
        target_trajs,
        bathymetry_data = bathy_data,
        strata = cfa_polys,
        title = "Snow Crab Larval Dispersal ($(opts.scenario), Year $(opts.projection_year))",
        output_path = fig1_path
    )

    # 2. DVM Depth Profiles
    fig2_path = joinpath(opts.output_dir, "dvm_depth_profiles.png")
    println("Rendering DVM depth profiles -> $(fig2_path)...")
    plot_vertical_migration_profiles(
        target_trajs,
        sample_indices = 1:min(5, size(target_trajs.depths, 1)),
        title = "Larval Diel Vertical Migration (DVM) Profiles",
        output_path = fig2_path
    )

    # 3. 2D Settlement Nursery Density Heatmap
    fig3_path = joinpath(opts.output_dir, "settlement_density.png")
    println("Rendering 2D settlement density distribution -> $(fig3_path)...")
    plot_larval_dispersal_density(
        target_trajs,
        title = "Snow Crab Nursery Settlement Density (%)",
        output_path = fig3_path
    )

    # 4. Empirical Movement Velocity Field and Diffusivity
    fig4_path = joinpath(opts.output_dir, "empirical_movement_field.png")
    println("Rendering empirical velocity quiver field -> $(fig4_path)...")
    emp_mov = estimate_empirical_movement(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_empirical_movement_field(emp_mov, output_path = fig4_path)

    # 5. Annotated Regional Connectivity Matrix
    fig5_path = joinpath(opts.output_dir, "regional_connectivity_matrix.png")
    println("Rendering regional connectivity matrix -> $(fig5_path)...")
    cfa_defs = if !isempty(cfa_polys)
        vcat(cfa_polys, [(name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))])
    else
        [
            (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
            (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
            (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
            (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
        ]
    end
    conn = compute_empirical_connectivity(target_trajs, strata_definitions = cfa_defs)
    plot_connectivity_matrix(conn, output_path = fig5_path)

    # 6. Thermal Exposure Map
    fig6_path = joinpath(opts.output_dir, "thermal_exposure_map.png")
    println("Rendering thermal exposure map -> $(fig6_path)...")
    therm_metrics = compute_gridded_thermal_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_thermal_exposure_map(therm_metrics, output_path = fig6_path)

    # 7. Recruitment Summary
    fig7_path = joinpath(opts.output_dir, "recruitment_summary.png")
    println("Rendering recruitment summary bar chart -> $(fig7_path)...")
    rec_metrics = compute_gridded_recruitment_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_recruitment_summary(rec_metrics, output_path = fig7_path)

    # 7b. Larval biology diagnostics: stage progression, thermal-time growth, and survival.
    # These were previously only embedded in the interactive HTML payload and never drawn, which
    # left the quantile/CDF stage model unverifiable from the run output.
    fig_molt_path = joinpath(opts.output_dir, "molt_progression.png")
    fig_dd_path = joinpath(opts.output_dir, "degree_day_growth.png")
    fig_surv_path = joinpath(opts.output_dir, "survival_curves.png")
    println("Rendering molt progression -> $(fig_molt_path)...")
    plot_molt_progression(target_trajs, output_path = fig_molt_path)
    println("Rendering degree-day growth -> $(fig_dd_path)...")
    plot_degree_day_growth(target_trajs, output_path = fig_dd_path)
    println("Rendering survival curves -> $(fig_surv_path)...")
    plot_survival_curves(target_trajs, output_path = fig_surv_path)

    # 8. Hydrodynamic Model Eulerian Advection & Tracers Figures
    fig8_path = joinpath(opts.output_dir, "hydrodynamic_advection.png")
    fig9_path = joinpath(opts.output_dir, "hydrodynamic_tracers.png")
    println("Rendering hydrodynamic advection velocity field -> $(fig8_path)...")
    active_hydro = if opts.enable_zarr && isdir(opts.zarr_path)
        try
            db_h = open_storage(opts.zarr_path; read_only = true)
            try
                target_run_id = !isempty(opts.run_id) ? opts.run_id :
                    "run_$(opts.scenario)_$(opts.projection_year)"
                load_hydrodynamic_field(db_h, target_run_id)
            finally
                close_storage(db_h)
            end
        catch
            nothing
        end
    else
        nothing
    end
    plot_hydrodynamic_advection(active_hydro, bathymetry_data = bathy_data, output_path = fig8_path)

    println("Rendering hydrodynamic seawater temperature & salinity tracers -> $(fig9_path)...")
    plot_hydrodynamic_tracers(active_hydro, output_path = fig9_path)

    fig_strat_path = joinpath(opts.output_dir, "hydrodynamic_stratification.png")
    println("Rendering hydrodynamic stratification diagnostics (N², S) -> $(fig_strat_path)...")
    plot_hydrodynamic_stratification(active_hydro, output_path = fig_strat_path)

    fig_diff_path = joinpath(opts.output_dir, "hydrodynamic_diffusion.png")
    println("Rendering hydrodynamic turbulent diffusion & eddy viscosity -> $(fig_diff_path)...")
    plot_hydrodynamic_diffusion(active_hydro, output_path = fig_diff_path)

    fig_sec_path = joinpath(opts.output_dir, "hydrodynamic_section.png")
    println("Rendering hydrodynamic vertical cross-section -> $(fig_sec_path)...")
    plot_hydrodynamic_section(
        active_hydro,
        variable = :temperature,
        coordinate = 44.0,
        section_type = :lat,
        output_path = fig_sec_path
    )

    # 8b. Eulerian time series at a representative station. This plot function existed but was
    # never called by the driver, so no run ever produced an Eulerian time series — the standard
    # check that a simulation is evolving plausibly at a fixed point.
    fig_ts_path = joinpath(opts.output_dir, "hydrodynamic_timeseries.png")
    if !isnothing(active_hydro)
        println("Rendering hydrodynamic station time series -> $(fig_ts_path)...")
        try
            plot_hydrodynamic_timeseries(
                active_hydro,
                station = (opts.domain_lon[1] + 0.55 * (opts.domain_lon[2] - opts.domain_lon[1]),
                           opts.domain_lat[1] + 0.55 * (opts.domain_lat[2] - opts.domain_lat[1])),
                variable = :temperature,
                output_path = fig_ts_path
            )
        catch err
            @warn "Station time-series plot failed (non-fatal): $(sprint(showerror, err))"
        end
    else
        fig_ts_path = nothing
    end

    # 8c. Additional Eulerian diagnostics. The time series above is now built from the real
    # snapshots; these three were absent entirely.
    fig_hov_path = joinpath(opts.output_dir, "hydrodynamic_hovmoller.png")
    fig_tbot_path = joinpath(opts.output_dir, "bottom_temperature_map.png")
    fig_ts_diagram_path = joinpath(opts.output_dir, "ts_diagram.png")
    if !isnothing(active_hydro)
        st_lon = opts.domain_lon[1] + 0.55 * (opts.domain_lon[2] - opts.domain_lon[1])
        st_lat = opts.domain_lat[1] + 0.55 * (opts.domain_lat[2] - opts.domain_lat[1])
        for (label, path, f) in (
            ("time-depth diagram", fig_hov_path,
             () -> plot_hydrodynamic_hovmoller(active_hydro; station = (st_lon, st_lat),
                                               variable = :temperature, output_path = fig_hov_path)),
            ("near-bed temperature map", fig_tbot_path,
             () -> plot_bottom_temperature_map(active_hydro; output_path = fig_tbot_path)),
            ("temperature-salinity diagram", fig_ts_diagram_path,
             () -> plot_temperature_salinity_diagram(active_hydro; output_path = fig_ts_diagram_path))
        )
            println("Rendering $label -> $(path)...")
            try
                f()
            catch err
                # A diagnostic that cannot be built from this archive must not abort the run.
                @warn "$label failed (non-fatal): $(first(split(sprint(showerror, err), "\n")))"
            end
        end
    else
        fig_hov_path = nothing
        fig_tbot_path = nothing
        fig_ts_diagram_path = nothing
    end

    # 9. Hydrodynamic Field & Dashboard Animation (Oceananigans / CairoMakie)
    anim_file_path = nothing
    if opts.animate_hydro
        anim_fmt = lowercase(opts.anim_format) == "gif" ? "gif" : "mp4"
        hydro_source = if !isnothing(active_hydro)
            active_hydro
        else
            default_jld2 = "hydrodynamics_$(opts.scenario)_$(opts.projection_year).jld2"
            resolved_jld2, _ = resolve_hydro_model_path(opts, default_jld2)
            isfile(resolved_jld2) ? resolved_jld2 : nothing
        end

        anim_file_path = if !isempty(strip(opts.anim_output_path))
            opts.anim_output_path
        elseif opts.anim_variable == :dashboard
            joinpath(opts.output_dir, "hydrodynamic_dashboard_animation.$(anim_fmt)")
        else
            joinpath(opts.output_dir, "hydrodynamic_$(opts.anim_variable)_animation.$(anim_fmt)")
        end

        if opts.anim_variable == :dashboard
            println("Rendering animated hydrodynamic dashboard -> $(anim_file_path)...")
            animate_hydrodynamic_dashboard(
                hydro_source;
                depth = opts.anim_depth,
                trajectories = target_trajs,
                bathymetry_data = bathy_data,
                output_path = anim_file_path,
                framerate = opts.anim_fps,
                show_trajectories = opts.anim_overlay_particles,
                domain_lon = opts.domain_lon,
                domain_lat = opts.domain_lat
            )
        else
            println("Rendering animated hydrodynamic $(opts.anim_variable) field -> $(anim_file_path)...")
            animate_hydrodynamic_field(
                hydro_source;
                variable = opts.anim_variable,
                depth = opts.anim_depth,
                trajectories = target_trajs,
                bathymetry_data = bathy_data,
                output_path = anim_file_path,
                framerate = opts.anim_fps,
                show_trajectories = opts.anim_overlay_particles,
                domain_lon = opts.domain_lon,
                domain_lat = opts.domain_lat
            )
        end
    end

    # 9. Query Archived Scenarios from Zarr for Cross-Scenario Comparison & Multi-Layer Interactive Map
    scenarios_bundle = Dict{String, Any}()
    if opts.enable_zarr && isdir(opts.zarr_path)
        try
            db = open_storage(opts.zarr_path; read_only = true)
            try
                runs_df = list_simulation_runs(db)
                for r in eachrow(runs_df)
                    s_id = string(r.run_id)
                    s_label = "$(r.scenario) ($(r.projection_year))"
                    try
                        s_trajs = load_trajectories_namedtuple(db, s_id)
                        s_disp = try load_gridded_dispersal(db, s_id) catch; nothing end
                        s_conn = try load_connectivity_matrix(db, s_id) catch; nothing end
                        s_hydro = try load_hydrodynamic_field(db, s_id) catch; nothing end
                        scenarios_bundle[s_label] = (
                            trajectories = s_trajs,
                            gridded_dispersal = s_disp,
                            connectivity = s_conn,
                            hydrodynamics = s_hydro
                        )
                    catch
                    end
                end
            finally
                close_storage(db)
            end
        catch
        end
    end

    if isempty(scenarios_bundle) && !isnothing(target_trajs)
        scenarios_bundle["$(opts.scenario) ($(opts.projection_year))"] = (
            trajectories = target_trajs,
            gridded_dispersal = (
                lon_centers = emp_mov.lon_centers,
                lat_centers = emp_mov.lat_centers,
                u_mean = emp_mov.u_mean,
                v_mean = emp_mov.v_mean,
                diffusivity = emp_mov.diffusivity,
                density = rec_metrics.settlement_density,
                mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                mean_degree_days = therm_metrics.mean_degree_days,
                sample_count = emp_mov.sample_count
            ),
            connectivity = conn,
            hydrodynamics = active_hydro
        )
    end

    # Cross-scenario 2D figure
    fig10_path = joinpath(opts.output_dir, "climate_scenario_comparison.png")
    println("Rendering cross-scenario climate comparison -> $(fig10_path)...")
    scenario_comp = Dict{Symbol, Any}()
    for (s_label, s_data) in scenarios_bundle
        scenario_comp[Symbol(s_label)] = s_data.trajectories
    end

    if length(scenario_comp) < 2
        flow_ssp585(lon, lat, z, t) = (0.09 + 0.03 * sin(t / 43200.0), -0.04, 0.0001)
        rng_comp = MersenneTwister(opts.seed + 100)
        larvae_comp = initialize_larval_particles(
            min(opts.n_particles, 50),
            lon_range = (opts.domain_lon[1] + 2.0, opts.domain_lon[2] - 3.0),
            lat_range = (opts.domain_lat[1] + 1.0, opts.domain_lat[2] - 1.0),
            min_seabed_depth = opts.min_seabed_depth,
            bathymetry = bathy_field_fn,
            rng = rng_comp
        )
        trajs_ssp585 = track_larval_cohort(
            larvae_comp,
            velocity_fn = flow_ssp585,
            total_duration = opts.track_duration,
            dt = opts.track_dt,
            # Same biology as the primary cohort, so the cross-scenario comparison is not
            # confounded by one branch running deterministic biology and the other stochastic.
            enable_molting = opts.enable_molting,
            cv_molt = opts.cv_molt,
            cv_mortality = opts.cv_mortality,
            cv_settlement = opts.cv_settlement,
            settlement_stochastic = opts.settlement_stochastic,
            rng = rng_comp
        )
        scenario_comp[:ssp585_2050] = trajs_ssp585
    end

    compare_scenario_dispersal(
        scenario_comp,
        title = "Scotian Shelf Snow Crab Dispersal Across Climate Scenarios",
        output_path = fig10_path
    )

    # 10. Standalone Multi-Layer Interactive HTML5 Dashboard
    html_path = joinpath(opts.output_dir, "interactive_larval_tracks.html")
    if opts.interactive_map
        println("Rendering multi-layer interactive Leaflet dashboard with hydrodynamic fields -> $(html_path)...")
        export_interactive_tracks_html(
            html_path;
            scenarios_data = scenarios_bundle,
            hydrodynamics = active_hydro,
            strata_definitions = cfa_defs,
            title = "Scotian Shelf Snow Crab Larval Dispersal & Demographic Connectivity"
        )
    end

    # 11. Depth-Stratified Voronoi Units Spatial Distribution
    fig_voronoi_path = nothing
    if opts.enable_voronoi
        fig_voronoi_path = joinpath(opts.output_dir, "voronoi_units_distribution.png")
        println("Rendering Voronoi units spatial distribution -> $(fig_voronoi_path)...")
        try
            # The run's own bathymetry, resolved the same way the grid was. The Voronoi units
            # are a partition of the same marine domain, so drawing them from a second product
            # would make the figure describe a different sea floor from the model it annotates.
            v_tess = generate_depth_stratified_voronoi_units(
                ;
                lon_range = opts.domain_lon,
                lat_range = opts.domain_lat,
                bathymetry_source = _run_bathymetry_source(opts),
                n_units = opts.voronoi_n_units,
                prob_core = opts.voronoi_prob_core,
                prob_shallow = opts.voronoi_prob_shallow,
                prob_deep = opts.voronoi_prob_deep,
                min_res_core_km = opts.voronoi_min_res_core_km,
                min_res_shallow_km = opts.voronoi_min_res_shallow_km,
                min_res_deep_km = opts.voronoi_min_res_deep_km,
                # TOML band edges are positive depths; the tessellator bins on signed z.
                core_depth_range = (-opts.voronoi_core_depth_max, -opts.voronoi_core_depth_min),
                shallow_depth_range = (-opts.voronoi_shallow_depth_max, -opts.voronoi_shallow_depth_min),
                deep_depth_range = (-opts.voronoi_deep_depth_max, -opts.voronoi_deep_depth_min),
                slope_weighting = opts.voronoi_slope_weighting,
                slope_factor = opts.voronoi_slope_factor,
                seed = opts.seed
            )
            fig_v = CairoMakie.Figure(size = (1000, 750), fontsize = 13)
            ax_v = CairoMakie.Axis(
                fig_v[1, 1],
                title = "Depth-Stratified Voronoi Units (N=$(length(v_tess.units)))",
                xlabel = "Longitude (°E)",
                ylabel = "Latitude (°N)"
            )
            lons = [u.lon for u in v_tess.units]
            lats = [u.lat for u in v_tess.units]
            strata_codes = [u.stratum == :core ? 1 : (u.stratum == :shallow ? 2 : 3) for u in v_tess.units]
            CairoMakie.scatter!(ax_v, lons, lats, color = strata_codes, colormap = :viridis, markersize = 6)
            CairoMakie.save(fig_voronoi_path, fig_v)
        catch v_err
            @warn "Failed to render Voronoi unit distribution figure: $(v_err)"
        end
    end

    println("All visualization figures and interactive maps successfully generated.")
    return (
        fig_trajectories = fig1_path,
        fig_dvm = fig2_path,
        fig_density = fig3_path,
        fig_empirical_movement = fig4_path,
        fig_connectivity = fig5_path,
        fig_thermal = fig6_path,
        fig_recruitment = fig7_path,
        # Larval biology diagnostics
        fig_molt_progression = fig_molt_path,
        fig_degree_day_growth = fig_dd_path,
        fig_survival = fig_surv_path,
        # Eulerian
        fig_hydro_advection = fig8_path,
        fig_hydro_tracers = fig9_path,
        fig_hydro_timeseries = fig_ts_path,
        fig_hovmoller = fig_hov_path,
        fig_bottom_temperature = fig_tbot_path,
        fig_ts_diagram = fig_ts_diagram_path,
        fig_comparison = fig10_path,
        fig_voronoi = fig_voronoi_path,
        interactive_map = opts.interactive_map ? html_path : nothing,
        animation = anim_file_path
    )
end


"""
    run_cli_list_runs(; opts::HydrodynamicOptions)

Query and display all simulation runs currently archived in Zarr.
"""
function run_cli_list_runs(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Archived Simulation Runs in Zarr: $(opts.zarr_path)")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    df = list_simulation_runs(db)
    close_storage(db)

    if nrow(df) == 0
        println("No simulation runs archived in database yet.")
        return
    end

    for row in eachrow(df)
        println("-----------------------------------------------------------------")
        println("Run ID:          $(row.run_id)")
        println("Scenario:        $(row.scenario) (Year $(row.projection_year))")
        println("Created:         $(row.created_at)")
        println("Particles:       $(row.n_particles) larvae | Duration: $(round(row.duration_days, digits=1)) days")
        println("Recruitment:     $(round(row.settlement_success_rate * 100, digits=2))% successful")
        println("Mean PLD:        $(round(row.mean_pld_days, digits=1)) days | Mean Temp: $(round(row.mean_exposure_temperature, digits=2)) °C")
        println("Dispersal Dist:  $(round(row.mean_dispersal_distance_km, digits=1)) km")
    end
    println("-----------------------------------------------------------------\n")
    return df
end


"""
    run_cli_compare_scenarios(; opts::HydrodynamicOptions)

Query and print comparative metrics across archived climate scenarios in Zarr.
"""
function run_cli_compare_scenarios(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Scenario Comparative Analysis from Zarr")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    df = compare_scenarios(db)
    close_storage(db)

    if nrow(df) == 0
        println("No simulation data available for scenario comparison.")
        return
    end

    for row in eachrow(df)
        println("Scenario: $(rpad(row.scenario, 16)) | Year: $(row.projection_year) | Runs: $(row.n_runs)")
        println("  Settlement Success: $(round(row.mean_settlement_success * 100, digits=2))%")
        println("  Mean PLD:           $(round(row.mean_pld_days, digits=1)) days")
        println("  Mean Exposure Temp: $(round(row.mean_temperature_celsius, digits=2)) °C")
        println("  Mean Dispersal:     $(round(row.mean_dispersal_km, digits=1)) km")
        println("  Thermal Mortality:  $(round(row.mean_thermal_mortality_rate * 100, digits=2))%")
        println()
    end
    return df
end


"""
    run_cli_model_average(; opts::HydrodynamicOptions)

Compute and display ensemble model-averaged demographic connectivity and recruitment.
"""
function run_cli_model_average(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Scenario Ensemble Model Averaging")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    runs_df = list_simulation_runs(db)

    if nrow(runs_df) == 0
        close_storage(db)
        println("No simulation data available for model averaging.")
        return
    end

    scens = unique(runs_df.scenario)
    println("Averaging across $(length(scens)) scenario models: $(join(scens, ", "))")
    ens = compute_ensemble_model_average(db, scens)
    close_storage(db)

    println("\nEnsemble Weighted Settlement Success: $(round(ens.mean_recruitment_rate * 100, digits=2))%")
    println("Ensemble Weighted Mean PLD:          $(round(ens.mean_pld_days, digits=1)) days")
    println("Ensemble Weighted Thermal Exposure:  $(round(ens.mean_thermal_exposure, digits=2)) °C")
    println("\nEnsemble Model-Averaged Connectivity Matrix (P_ij ± std):")
    n_s = length(ens.strata_names)
    for i in 1:n_s
        row_strs = [
            "$(round(ens.mean_connectivity[i, j], digits=3))±$(round(ens.std_connectivity[i, j], digits=3))"
            for j in 1:n_s
        ]
        println("  $(rpad(ens.strata_names[i], 28)) -> [$(join(row_strs, ", "))]")
    end
    println()
    return ens
end

"""
    animate_hydrodynamic_dashboard(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
animate_hydrodynamic_dashboard(args...; kwargs...) = error("`animate_hydrodynamic_dashboard` is owned by GeoOceanModel, not by this package.")

"""
    animate_hydrodynamic_field(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
animate_hydrodynamic_field(args...; kwargs...) = error("`animate_hydrodynamic_field` is owned by GeoOceanModel, not by this package.")

"""
    plot_hydrodynamic_advection(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
plot_hydrodynamic_advection(args...; kwargs...) = error("`plot_hydrodynamic_advection` is owned by GeoOceanModel, not by this package.")

"""
    plot_hydrodynamic_diffusion(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
plot_hydrodynamic_diffusion(args...; kwargs...) = error("`plot_hydrodynamic_diffusion` is owned by GeoOceanModel, not by this package.")

"""
    plot_hydrodynamic_section(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
plot_hydrodynamic_section(args...; kwargs...) = error("`plot_hydrodynamic_section` is owned by GeoOceanModel, not by this package.")

"""
    plot_hydrodynamic_stratification(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
plot_hydrodynamic_stratification(args...; kwargs...) = error("`plot_hydrodynamic_stratification` is owned by GeoOceanModel, not by this package.")

"""
    plot_hydrodynamic_timeseries(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
plot_hydrodynamic_timeseries(args...; kwargs...) = error("`plot_hydrodynamic_timeseries` is owned by GeoOceanModel, not by this package.")

"""
    plot_hydrodynamic_tracers(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
plot_hydrodynamic_tracers(args...; kwargs...) = error("`plot_hydrodynamic_tracers` is owned by GeoOceanModel, not by this package.")

"""
    _run_bathymetry(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
_run_bathymetry(args...; kwargs...) = error("`_run_bathymetry` is owned by GeoOceanModel, not by this package.")

"""
    _run_bathymetry_source(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
_run_bathymetry_source(args...; kwargs...) = error("`_run_bathymetry_source` is owned by GeoOceanModel, not by this package.")

"""
    expand_domain_with_buffer(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
expand_domain_with_buffer(args...; kwargs...) = error("`expand_domain_with_buffer` is owned by GeoOceanModel, not by this package.")

"""
    flow_ssp585(args...; kwargs...)

Owned by the OCEAN half. GeoLarvalTracking calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
flow_ssp585(args...; kwargs...) = error("`flow_ssp585` is owned by GeoOceanModel, not by this package.")
