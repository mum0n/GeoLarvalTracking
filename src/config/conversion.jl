"""
    conversion.jl

Configuration conversion.

Turning the nested dictionary into runtime options and into the three pure config
structs, and back. `get_val` lives here because it is how a dict may be missing
keys during conversion.
"""

"""
    _PROJECT_ROOT

The root of this package -- the directory holding `src/`, `configs/` and `Project.toml`.

Configured output paths are resolved against it, so a run started from any working
directory writes to the same place.
"""
const _PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))


"""
    configuration_to_options(config_dict::AbstractDict; config_name = "", overrides...) -> HydrodynamicOptions

Convert a configuration dictionary into a validated `HydrodynamicOptions` instance,
allowing optional keyword parameter overrides.

# Inputs
- `config_dict::AbstractDict`: Parsed configuration dictionary.
- `config_name::AbstractString`: The configuration's own name, used to name its directory
  in the lakehouse when the configuration names no output directory of its own. Empty
  means "the caller did not say", and the name `default` is used.
- `overrides...`: Additional keyword arguments to override dictionary entries.

# Outputs
- `HydrodynamicOptions`: Constructed runtime options instance.
"""
function configuration_to_options(config_dict::AbstractDict;
                                  config_name::AbstractString = "",
                                  overrides...)
    # Helper to safely extract nested values with defaults
    function get_val(section::String, key::String, default_val)
        if haskey(config_dict, section) && haskey(config_dict[section], key)
            return config_dict[section][key]
        end
        return default_val
    end

    lon_min = Float64(get_val("domain", "lon_min", -71.0))
    lon_max = Float64(get_val("domain", "lon_max", -53.0))
    lat_min = Float64(get_val("domain", "lat_min", 40.0))
    lat_max = Float64(get_val("domain", "lat_max", 48.5))
    z_min   = Float64(get_val("domain", "z_min", -3500.0))
    z_max   = Float64(get_val("domain", "z_max", 0.0))

    buffer_km = Float64(get_val("domain", "buffer_km", get_val("biology", "buffer_km", 100.0)))

    nx = Int(get_val("grid", "nx", 50))
    ny = Int(get_val("grid", "ny", 50))
    nz = Int(get_val("grid", "nz", 10))

    bathy_source = String(get_val("data", "bathy_source", "etopo2022"))
    coastline_source = String(get_val("data", "coastline_source", "natural_earth_10m"))
    wind_source = String(get_val("data", "wind_source", "open_meteo"))
    wind_time_iso = String(get_val("data", "wind_time_iso", "2020-06-01T00:00:00Z"))

    # Shelf geometry, from [bathymetry]. `inshore_depth` is the water depth at the landward
    # edge of the domain; `shelf_slope` is the total rise across the shelf. Both are read
    # from here rather than hardcoded, so a config can actually shape its own seabed.
    inshore_depth = Float64(get_val("bathymetry", "inshore_depth", -20.0))
    shelf_slope   = Float64(get_val("bathymetry", "shelf_slope", 500.0))
    mask_bay_of_fundy = Bool(get_val("bathymetry", "mask_bay_of_fundy", false))

    # an empty read is not an error; it simply means no constituents were requested.
    #
    # case-sensitively, while the atlas file lookup is case-insensitive. Normalising here


    sim_hours = Float64(get_val("hydrodynamics", "sim_duration_hours", 12.0))
    sim_dt    = Float64(get_val("hydrodynamics", "sim_dt_seconds", 120.0))
    # Floor on the adaptive time step, forwarded to the simulation as `min_Δt`. The surface cell
    # sets the admissible step (Δt ≈ target_cfl × Δz_min / w), so a floor above that value stops
    # the wizard stabilising the surface layer.
    min_dt_seconds = Float64(get_val("hydrodynamics", "min_dt_seconds", 2.0))
    adapt_cfl = Bool(get_val("hydrodynamics", "adaptive_cfl", true))
    tgt_cfl   = Float64(get_val("hydrodynamics", "target_cfl", 0.2))
    heat_flux = Float64(get_val("hydrodynamics", "surface_heat_flux", 50.0))
    # `bulk_heat_flux` selects WHICH heat flux is applied: true computes it from the ingested
    # atmospheric surface state, false applies the constant `surface_heat_flux` above. It was
    # parsed into `HydrodynamicConfig` and documented in the TOML but never carried onto
    # `HydrodynamicOptions`, so nothing downstream could read it and the bulk flux was applied
    # unconditionally.
    #
    # It is read from `[atmosphere]`, which is where the TOML actually carries it (next to
    # `atmospheric_source`, `wind_source` and the rest of the atmospheric configuration). The
    # `[hydrodynamics]` section is accepted as a fallback so the key is not silently ignored
    # if it is written there instead.
    bulk_heat_flux = if haskey(config_dict, "atmosphere") &&
                        haskey(config_dict["atmosphere"], "bulk_heat_flux")
        Bool(get_val("atmosphere", "bulk_heat_flux", true))
    else
        Bool(get_val("hydrodynamics", "bulk_heat_flux", true))
    end
    # Cap on how many hydrodynamic snapshots are held in memory when building the 4D flow
    # interpolator for particle tracking. A production run stores thousands of snapshots and
    # materialising all of them for u, v, w and T is hundreds of gigabytes; the snapshots are
    # evenly subsampled down to this count. `0` means "no cap" (load every snapshot).
    max_flow_snapshots = Int(get_val("hydrodynamics", "max_flow_snapshots", 400))
    allow_analytical_fallback = Bool(get_val("hydrodynamics", "allow_analytical_fallback", false))
    # Cap on the background current a larva is advected by. A diverged hydrodynamics cell can
    # otherwise teleport a particle across the domain in a single step; 3 m/s is well above any
    # real Scotian Shelf current. Use a negative value to disable the clamp.
    _mcs = Float64(get_val("hydrodynamics", "max_current_speed", 3.0))
    max_current_speed = _mcs > 0 ? _mcs : nothing

    n_parts   = Int(get_val("biology", "n_particles", 100))
    track_days = Float64(get_val("biology", "track_duration_days", 5.0))
    track_dt  = Float64(get_val("biology", "track_dt_seconds", 300.0))
    min_depth = Float64(get_val("biology", "min_seabed_depth", 100.0))
    diff_h    = Float64(get_val("biology", "diffusivity_h", 10.0))
    diff_v    = Float64(get_val("biology", "diffusivity_v", 1e-4))

    rel_mode_str = String(get_val("biology", "release_depth_mode", "bottom"))
    rel_mode = Symbol(lowercase(rel_mode_str))
    raw_offset = get_val("biology", "bottom_release_offset", [0.5, 3.0])
    b_offset = if raw_offset isa Vector && length(raw_offset) >= 2
        (Float64(raw_offset[1]), Float64(raw_offset[2]))
    elseif raw_offset isa Tuple && length(raw_offset) >= 2
        (Float64(raw_offset[1]), Float64(raw_offset[2]))
    else
        (0.5, 3.0)
    end
    init_ascent = Bool(get_val("biology", "enable_initial_ascent", true))
    asc_spd     = Float64(get_val("biology", "ascent_speed", 0.010))
    asc_target  = Float64(get_val("biology", "ascent_target_depth", -10.0))

    enable_dvm = Bool(get_val("dvm", "enable_dvm", true))
    enable_molting = Bool(get_val("molting_and_settlement", "enable_molting", true))

    # Diel Vertical Migration: the target depths the larvae hold by day and by night,
    # per stage, and the two swim-speed multipliers. These decide where the cohort is in
    # the water column at any hour, so they are read from the configuration rather than
    # left to the behaviour module's defaults -- a run whose DVM depths live in code
    # cannot be reproduced from its configuration file.
    zoea1_day_depth       = Float64(get_val("dvm", "zoea1_day_depth", -50.0))
    zoea1_night_depth     = Float64(get_val("dvm", "zoea1_night_depth", -10.0))
    zoea2_day_depth       = Float64(get_val("dvm", "zoea2_day_depth", -55.0))
    zoea2_night_depth     = Float64(get_val("dvm", "zoea2_night_depth", -8.0))
    megalopa_day_depth    = Float64(get_val("dvm", "megalopa_day_depth", -120.0))
    megalopa_night_depth  = Float64(get_val("dvm", "megalopa_night_depth", -60.0))
    zoea2_depth_factor    = Float64(get_val("dvm", "zoea2_depth_factor", 1.2))
    megalopa_swim_factor  = Float64(get_val("dvm", "megalopa_swim_factor", 1.5))

    # Phenology: the degree-day thresholds that move a larva between stages, and the
    # mortality and settlement envelopes. A cohort that molts at 65 °C·d in one run and
    # 150 °C·d in the next is a different model, and the difference has to be in the
    # configuration file for anyone to see it.
    t_base                   = Float64(get_val("molting_and_settlement", "t_base", -1.5))
    dd_zoea1_to_zoea2        = Float64(get_val("molting_and_settlement", "dd_zoea1_to_zoea2", 65.0))
    dd_zoea2_to_megalopa     = Float64(get_val("molting_and_settlement", "dd_zoea2_to_megalopa", 130.0))
    dd_megalopa_to_settle    = Float64(get_val("molting_and_settlement", "dd_megalopa_to_settle", 200.0))
    mortality_base           = Float64(get_val("molting_and_settlement", "mortality_base", 0.02))
    mortality_thermal_threshold = Float64(get_val("molting_and_settlement", "mortality_thermal_threshold", 7.0))
    mortality_thermal_sensitivity = Float64(get_val("molting_and_settlement", "mortality_thermal_sensitivity", 0.015))
    mortality_cold_threshold = Float64(get_val("molting_and_settlement", "mortality_cold_threshold", -1.5))
    mortality_cold_sensitivity = Float64(get_val("molting_and_settlement", "mortality_cold_sensitivity", 0.01))
    settlement_min_depth     = Float64(get_val("molting_and_settlement", "settlement_min_depth", -250.0))
    settlement_max_depth     = Float64(get_val("molting_and_settlement", "settlement_max_depth", -50.0))
    settlement_max_temp      = Float64(get_val("molting_and_settlement", "settlement_max_temp", 6.0))

    # Stochastic biology parameters
    settlement_stochastic = Bool(get_val("biology", "settlement_stochastic", true))
    cv_molt = Float64(get_val("biology", "cv_molt", 0.25))
    cv_mortality = Float64(get_val("biology", "cv_mortality", 0.25))
    cv_settlement = Float64(get_val("biology", "cv_settlement", 0.25))

    interactive  = Bool(get_val("visualization", "interactive_map", true))
    anim_hydro   = Bool(get_val("visualization", "animate_hydro", false))
    anim_var     = Symbol(lowercase(String(get_val("visualization", "anim_variable", "dashboard"))))
    anim_fps     = Int(get_val("visualization", "anim_fps", 10))
    anim_fmt     = String(get_val("visualization", "anim_format", "mp4"))
    anim_depth   = Float64(get_val("visualization", "anim_depth", -2.5))
    anim_overlay = Bool(get_val("visualization", "anim_overlay_particles", false))
     anim_out_path = String(get_val("visualization", "anim_output_path", ""))

     # Where this half writes, and what it reads. Both are storage: the directory the run
     # writes into, and the realized ocean store it tracks through. One section, one
     # source of truth -- two sections naming the same run's places is two chances for
     # them to disagree, and the disagreement is silent: the run writes its results
     # somewhere the next run does not look.
     _storage_output = String(get_val("storage", "output_dir", ""))
     # Where a run writes when the configuration names no directory of its own: the
     # lakehouse, under this configuration's own name. NOT a package-relative
     # `outputs/default`, which put the run in the package root -- outside the lakehouse,
     # where the other half never looks and no later run finds it. `lakehouse_dir` is the
     # shared resolver, so the environment override, the local preference and the user
     # configuration apply here as they do everywhere else.
     _config_stem = isempty(strip(config_name)) ? "default" : strip(String(config_name))
     _lakehouse_output = GeoData.lakehouse_dir("GeoLarvalTracking", _config_stem)
     # The path is resolved against THIS PACKAGE'S ROOT, not the current working directory:
     # a run started from anywhere then lands in the same place. The resolver also accepts
     # the lakehouse either inside the project (`lakehouse/...`) or beside it
     # (`../lakehouse/...`), so a config written for one layout still resolves in the other.
     # When the configuration names no directory at all, the lakehouse default is used --
     # resolved absolutely, so a run started from any working directory writes to the
     # same place.
     output_dir = isempty(_storage_output) ? _lakehouse_output :
                  GeoData.project_lakehouse_path(_storage_output, _PROJECT_ROOT)
     # The trajectory store lives under this run's output directory unless the
     # configuration names a path explicitly. A CWD-relative path put the store
     # outside the lakehouse whenever the run was started from elsewhere.
     _zarr_raw = String(get_val("storage", "zarr_path", ""))
     storage_backend = String(get_val("storage", "storage_backend", "zarr"))
     enable_zarr = Bool(get_val("storage", "enable_zarr", true))
     
     enable_checkpoint = Bool(get_val("storage", "enable_checkpoint", true))
     checkpoint_prefix = String(get_val("storage", "checkpoint_prefix", "checkpoint_ParticleTracking"))
     checkpoint_schedule_seconds = Float64(get_val("storage", "checkpoint_schedule_seconds", 21600.0))
     checkpoint_dir = String(get_val("storage", "checkpoint_dir", ""))
     checkpoint_cleanup = Bool(get_val("storage", "checkpoint_cleanup", true))
     
     # Ingested data lives in the run's output directory's `inputs` subdirectory, beside
     # this run's outputs, provenance and resolved config. It is deliberately NOT
     # configurable and NOT shared between runs: a single top-level `inputs/` let a run
     # with a changed domain silently reuse another domain's bathymetry, which is exactly
     # the quiet mismatch `resolved_config.toml` and `data_provenance.json` exist to
     # expose. The cost is disk -- a dataset is fetched once per configuration rather than once
     # per machine -- and that is the intended trade, because it makes each configuration a
     # self-contained, shareable directory.
     input_dir  = joinpath(output_dir, "inputs")
     # The trajectory store and the restart archives are part of this run's output
     # tree. Named explicitly, they are honoured; otherwise they are derived.
     zarr_path = isempty(strip(_zarr_raw)) ?
                 joinpath(output_dir, "particle_tracking.zarr") :
                 GeoData.project_lakehouse_path(_zarr_raw, _PROJECT_ROOT)
     seed       = Int(get_val("runtype", "seed", 42))

     # The realized ocean store this run tracks through. `[storage] ocean_model` names it;
     # `[hydrodynamics] hydro_model_file` is accepted as an alias. When neither names one,
     # the ocean half's own lakehouse directory is the default -- the ocean run of the
     # same name, not a file inside this run's output directory, which is where an
     # unqualified fallback used to look and where no ocean store has ever been.
     hydro_file = String(get_val("storage", "ocean_model",
                          get_val("hydrodynamics", "hydro_model_file",
                              GeoData.lakehouse_dir("GeoOceanModel", _config_stem))))
     hydro_only = Bool(get_val("hydrodynamics", "hydro_only", false))
     track_only = Bool(get_val("hydrodynamics", "track_only", false))
    reuse_hydro = Bool(get_val("hydrodynamics", "reuse_hydro", false))
    run_id_val = String(get_val("storage", "run_id", ""))

    enable_cp = Bool(get_val("storage", "enable_checkpoint",
                             get_val("hydrodynamics", "enable_checkpoint", true)))
    cp_pfx_raw = String(get_val("storage", "checkpoint_prefix",
                                get_val("hydrodynamics", "checkpoint_prefix", "")))
    cp_sched  = Float64(get_val("storage", "checkpoint_schedule_seconds",
                                get_val("hydrodynamics", "checkpoint_schedule_seconds", 0.0)))
    # Field-output cadence, distinct from the checkpoint cadence. Read from [storage] rather
    # than hardcoded downstream, so the two can legitimately differ.
    out_sched = Float64(get_val("storage", "output_schedule_seconds", 21600.0))
    # Restart archives are part of this run's output tree. Named explicitly they
    # are honoured; otherwise they sit beside the trajectory store.
    _cp_raw = String(get_val("storage", "checkpoint_dir", ""))
    cp_dir = isempty(strip(_cp_raw)) ?
               joinpath(output_dir, "checkpoints") :
               GeoData.project_lakehouse_path(_cp_raw, _PROJECT_ROOT)
    cp_clean  = Bool(get_val("storage", "checkpoint_cleanup", true))
    auto_res  = Bool(get_val("hydrodynamics", "auto_restart", true))
    overwrite_val = Bool(get_val("hydrodynamics", "overwrite", false))
    v_mode_str = String(get_val("grid", "vertical_stretching_mode", "oceananigans"))
    # Two-segment shelf grid parameters. `vertical_stretching_mode = "two_segment"` splits `nz`
    # into a surface-refined upper segment above `vertical_break_depth` and a coarse lower segment,
    # which is the only way to resolve the 10-400 m active layer while still reaching the deep basin.
    # `nz_above = 0` means "derive it as 60% of nz".
    nz_above_cfg = Int(get_val("grid", "nz_above", 0))
    vertical_break_depth = Float64(get_val("grid", "vertical_break_depth", -400.0))
    v_mode = Symbol(lowercase(v_mode_str))
    v_file = strip(String(get_val("grid", "vertical_grid_file", "")))

    # `vertical_grid_file` is the operator's own file: a relative path is resolved against the
    # working directory they ran from, not against a previous run's inputs. If they name a
    # file, it has to exist -- silently falling back to a generated grid would run the whole
    # simulation on a vertical grid nobody chose. The file is then copied into this run's
    # containerised input directory so the run is self-contained and reproducible.
    #
    # `[grid] depths` is the discretised vertical axis written inline in the config, as
    # positive-downward face depths starting at the surface: depths = [0, 10, 20, ...].
    # It is the alternative to naming a CSV, and the two are mutually exclusive so a config
    # can never silently prefer one over the other.
    raw_depths = get_val("grid", "depths", Any[])
    v_depths = Float64[Float64(d) for d in raw_depths]
    if !isempty(v_depths)
        isempty(v_file) || error(
            "[grid] sets both vertical_grid_file and depths. Use one: an explicit depth list " *
            "is easier to review in the config than a side-car CSV.")
        abs(first(v_depths)) <= 1e-6 || error(
            "[grid] depths must start at 0.0 (the surface); got $(first(v_depths)).")
        all(diff(v_depths) .> 0) || error(
            "[grid] depths must increase monotonically (positive downward): $(v_depths).")
        length(v_depths) == nz + 1 || error(
            "[grid] depths has $(length(v_depths)) entries, but nz = $nz requires " *
            "$(nz + 1) faces.")
        # The deepest face has to reach the bottom of the domain, or the grid silently stops
        # short of the seabed and the deepest water column is left unresolvable.
        abs(last(v_depths) - abs(z_min)) <= 1.0 || error(
            "[grid] depths ends at $(last(v_depths)) m but [domain] z_min is $(z_min) m. " *
            "The deepest face must reach the domain bottom.")
    end
    if !isempty(v_file)
        src = abspath(v_file)
        isfile(src) || error(
            "[grid] vertical_grid_file names '$v_file', which does not exist (resolved to " *
            "$src). Supply the file, or remove the key to derive the grid from the depth " *
            "range and nz instead."
        )
        dst = joinpath(input_dir, basename(v_file))
        isfile(dst) || cp(src, dst, force = true)
        v_file = dst
    end
    atmo_src = Symbol(lowercase(String(get_val("atmosphere", "forcing", "era5_climatology"))))
    obc_src_val = get_val("boundaries", "ocean_boundary_source", nothing)
    if isnothing(obc_src_val)
        obc_src_val = get_val("boundaries", "parent_ocean", "glorys12v1")
    end
    obc_src = Symbol(lowercase(String(obc_src_val)))
    obc_tp = Symbol(lowercase(String(get_val("boundaries", "obc_type", "flather_chapman"))))

    # Initial 3-D hydrography (T, S). Read from `[data]`, independent of the
    # bathymetry/wind `bathy_source` and the lateral `ocean_boundary_source`, because the
    # initial water column is a separate choice from the forcing and the boundaries.
    hydro_src = Symbol(lowercase(String(get_val("data", "hydrography_source", "synthetic"))))
    hydro_month = Int(get_val("data", "hydrography_month", 0))

    # Lateral boundary relaxation (Price & Aumont 2011) parameters
    boundary_method = Symbol(lowercase(String(get_val("boundaries", "method", "relaxation"))))
    sponge_width_val = Float64(get_val("boundaries", "sponge_layer_width", 0.35))
    sponge_tau_val   = Float64(get_val("boundaries", "sponge_timescale_seconds", 3600.0))
    u_inflow_val     = Float64(get_val("boundaries", "u_inflow", -0.15))
    v_inflow_val     = Float64(get_val("boundaries", "v_inflow", 0.05))

    # Maximum allowable time step ceiling
    max_dt_val = Float64(get_val("hydrodynamics", "max_dt_seconds", sim_dt))

    # Construct HydrodynamicOptions with overrides applied
    return HydrodynamicOptions(;
        domain_lon = (lon_min, lon_max),
        domain_lat = (lat_min, lat_max),
        domain_z = (z_min, z_max),
        grid_size = (nx, ny, nz),
        bathy_source = bathy_source,
    coastline_source = coastline_source,
    wind_source = wind_source,
    wind_time_iso = wind_time_iso,
        inshore_depth = inshore_depth,
        shelf_slope = shelf_slope,
        sim_dt = sim_dt,
        sim_duration = sim_hours * 3600.0,
        min_dt_seconds = min_dt_seconds,
        adaptive_cfl = adapt_cfl,
        target_cfl = tgt_cfl,
          surface_heat_flux = heat_flux,
          bulk_heat_flux = Bool(bulk_heat_flux),
    max_flow_snapshots = max_flow_snapshots,
    allow_analytical_fallback = allow_analytical_fallback,
    max_current_speed = (max_current_speed isa Real && max_current_speed > 0) ?
                        Float64(max_current_speed) : nothing,
        n_particles = n_parts,
        track_duration = track_days * 86400.0,
        track_dt = track_dt,
        diffusivity_h = diff_h,
        diffusivity_v = diff_v,
        enable_dvm = enable_dvm,
        enable_molting = enable_molting,
        zoea1_day_depth = zoea1_day_depth,
        zoea1_night_depth = zoea1_night_depth,
        zoea2_day_depth = zoea2_day_depth,
        zoea2_night_depth = zoea2_night_depth,
        megalopa_day_depth = megalopa_day_depth,
        megalopa_night_depth = megalopa_night_depth,
        zoea2_depth_factor = zoea2_depth_factor,
        megalopa_swim_factor = megalopa_swim_factor,
        t_base = t_base,
        dd_zoea1_to_zoea2 = dd_zoea1_to_zoea2,
        dd_zoea2_to_megalopa = dd_zoea2_to_megalopa,
        dd_megalopa_to_settle = dd_megalopa_to_settle,
        mortality_base = mortality_base,
        mortality_thermal_threshold = mortality_thermal_threshold,
        mortality_thermal_sensitivity = mortality_thermal_sensitivity,
        mortality_cold_threshold = mortality_cold_threshold,
        mortality_cold_sensitivity = mortality_cold_sensitivity,
        settlement_min_depth = settlement_min_depth,
        settlement_max_depth = settlement_max_depth,
        settlement_max_temp = settlement_max_temp,
        min_seabed_depth = min_depth,
        buffer_km = buffer_km,
        release_depth_mode = rel_mode,
        bottom_release_offset = b_offset,
        enable_initial_ascent = init_ascent,
        ascent_speed = asc_spd,
        ascent_target_depth = asc_target,
        use_gpu = false,
        fallback_to_cpu = false,
          interactive_map = interactive,
        zarr_path = zarr_path,
        storage_backend = storage_backend,
        enable_zarr = enable_zarr,
        output_dir = output_dir,
        input_dir = input_dir,
        seed = seed,
        hydro_model_file = hydro_file,
        hydro_only = hydro_only,
        track_only = track_only,
        reuse_hydro = reuse_hydro,
        enable_checkpoint = enable_cp,
        checkpoint_prefix = cp_pfx_raw,
        checkpoint_schedule = cp_sched,
        output_schedule_seconds = out_sched,
        checkpoint_dir = cp_dir,
        checkpoint_cleanup = cp_clean,
        auto_restart = auto_res,
        overwrite = overwrite_val,
        run_id = run_id_val,
        vertical_stretching_mode = v_mode,
        nz_above = nz_above_cfg,
        vertical_break_depth = vertical_break_depth,
        vertical_grid_file = v_file,
        vertical_depths = v_depths,
        atmospheric_source = atmo_src,
    ocean_boundary_source = obc_src,
    embedding_lon = embedding_domain_ranges(config_dict).lon,
    embedding_lat = embedding_domain_ranges(config_dict).lat,
    hydrography_source = hydro_src,
    hydrography_month = hydro_month,
        obc_type = obc_tp,
        animate_hydro = anim_hydro,
        anim_variable = anim_var,
        anim_fps = anim_fps,
        anim_format = anim_fmt,
        anim_depth = anim_depth,
        anim_overlay_particles = anim_overlay,
        anim_output_path = anim_out_path,
        boundary_method = boundary_method,
        sponge_layer_width = sponge_width_val,
        sponge_timescale = sponge_tau_val,
        u_inflow = u_inflow_val,
        v_inflow = v_inflow_val,
        max_dt = max_dt_val,
        settlement_stochastic = settlement_stochastic,
        cv_molt = cv_molt,
        cv_mortality = cv_mortality,
        cv_settlement = cv_settlement,
        mask_bay_of_fundy = mask_bay_of_fundy,
        overrides...
    )
end

"""
    options_to_configuration(opts::HydrodynamicOptions) -> Dict{String, Any}

Extract a complete configuration dictionary from a `HydrodynamicOptions` instance.
"""
function options_to_configuration(opts::HydrodynamicOptions)::Dict{String, Any}
    return Dict{String, Any}(
        "domain" => Dict{String, Any}(
            "lon_min" => opts.domain_lon[1],
            "lon_max" => opts.domain_lon[2],
            "lat_min" => opts.domain_lat[1],
            "lat_max" => opts.domain_lat[2],
            "z_min" => opts.domain_z[1],
            "z_max" => opts.domain_z[2],
            "buffer_km" => opts.buffer_km
        ),
        "grid" => Dict{String, Any}(
            "nx" => opts.grid_size[1],
            "ny" => opts.grid_size[2],
            "nz" => opts.grid_size[3],
            "vertical_stretching_mode" => string(opts.vertical_stretching_mode),
            "nz_above" => opts.nz_above,
            "vertical_break_depth" => opts.vertical_break_depth,
            "vertical_grid_file" => opts.vertical_grid_file,
            "depths" => opts.vertical_depths
        ),
        "atmosphere" => Dict{String, Any}(
            "source" => string(opts.atmospheric_source)
        ),
        "boundaries" => Dict{String, Any}(
            "ocean_boundary_source" => string(opts.ocean_boundary_source),
            "obc_type" => string(opts.obc_type),
            "method" => string(opts.boundary_method),
            "sponge_layer_width" => opts.sponge_layer_width,
            "sponge_timescale_seconds" => opts.sponge_timescale,
            "u_inflow" => opts.u_inflow,
            "v_inflow" => opts.v_inflow
        ),
        "data" => Dict{String, Any}(
            "bathy_source" => opts.bathy_source,
            "coastline_source" => opts.coastline_source,
            "wind_source" => opts.wind_source,
            "wind_time_iso" => opts.wind_time_iso,
            "hydrography_source" => string(opts.hydrography_source),
            "hydrography_month" => opts.hydrography_month
        ),
        "bathymetry" => Dict{String, Any}(
            "inshore_depth" => opts.inshore_depth,
            "shelf_slope" => opts.shelf_slope,
            "mask_bay_of_fundy" => opts.mask_bay_of_fundy
        ),
        "tides" => Dict{String, Any}(
        ),
        "hydrodynamics" => Dict{String, Any}(
            "sim_duration_hours" => opts.sim_duration / 3600.0,
            "sim_dt_seconds" => opts.sim_dt,
            "adaptive_cfl" => opts.adaptive_cfl,
            "target_cfl" => opts.target_cfl,
            "max_dt_seconds" => opts.max_dt,
            "min_dt_seconds" => opts.min_dt_seconds,
            "surface_heat_flux" => opts.surface_heat_flux,
            "max_flow_snapshots" => something(opts.max_flow_snapshots, 0),
            "allow_analytical_fallback" => opts.allow_analytical_fallback,
            "max_current_speed" => something(opts.max_current_speed, 0.0),
            "hydro_model_file" => opts.hydro_model_file,
            "hydro_only" => opts.hydro_only,
            "track_only" => opts.track_only,
            "reuse_hydro" => opts.reuse_hydro,
            "overwrite" => opts.overwrite
        ),
        "biology" => Dict{String, Any}(
            "n_particles" => opts.n_particles,
            "track_duration_days" => opts.track_duration / 86400.0,
            "track_dt_seconds" => opts.track_dt,
            "min_seabed_depth" => opts.min_seabed_depth,
            "buffer_km" => opts.buffer_km,
            "diffusivity_h" => opts.diffusivity_h,
            "diffusivity_v" => opts.diffusivity_v,
            "release_depth_mode" => string(opts.release_depth_mode),
            "bottom_release_offset" => [opts.bottom_release_offset[1], opts.bottom_release_offset[2]],
            "enable_initial_ascent" => opts.enable_initial_ascent,
            "ascent_speed" => opts.ascent_speed,
            "ascent_target_depth" => opts.ascent_target_depth,
            "settlement_stochastic" => opts.settlement_stochastic,
            "cv_molt" => opts.cv_molt,
            "cv_mortality" => opts.cv_mortality,
            "cv_settlement" => opts.cv_settlement
        ),
        "dvm" => Dict{String, Any}(
            "enable_dvm" => opts.enable_dvm
        ),
        "molting_and_settlement" => Dict{String, Any}(
            "enable_molting" => opts.enable_molting,
            "cv_molt" => opts.cv_molt,
            "cv_mortality" => opts.cv_mortality
        ),
        "storage" => Dict{String, Any}(
            "output_dir" => opts.output_dir,
            "ocean_model" => opts.hydro_model_file,
            "enable_zarr" => opts.enable_zarr,
            "zarr_path" => opts.zarr_path,
            "run_id" => opts.run_id,
            "enable_checkpoint" => opts.enable_checkpoint,
            "checkpoint_prefix" => opts.checkpoint_prefix,
            "checkpoint_schedule_seconds" => opts.checkpoint_schedule,
            "output_schedule_seconds" => opts.output_schedule_seconds,
            "checkpoint_dir" => opts.checkpoint_dir,
            "checkpoint_cleanup" => opts.checkpoint_cleanup
        ),
        "runtype" => Dict{String, Any}(
            "seed" => opts.seed
        ),
        "visualization" => Dict{String, Any}(
            "interactive_map" => opts.interactive_map,
            "animate_hydro" => opts.animate_hydro,
            "anim_variable" => string(opts.anim_variable),
            "anim_fps" => opts.anim_fps,
            "anim_format" => opts.anim_format,
            "anim_depth" => opts.anim_depth,
            "anim_overlay_particles" => opts.anim_overlay_particles,
            "anim_output_path" => opts.anim_output_path
        )
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Formally Decoupled Configuration Architecture
# ─────────────────────────────────────────────────────────────────────────────


"""
    to_hydrodynamic_config(opts::HydrodynamicOptions; kwargs...) -> HydrodynamicConfig

Extract and convert a `HydrodynamicOptions` instance into a pure `HydrodynamicConfig`.
"""
function to_hydrodynamic_config(opts::HydrodynamicOptions; kwargs...)::HydrodynamicConfig
    d = Dict{Symbol, Any}(
        :domain_lon               => opts.domain_lon,
        :domain_lat               => opts.domain_lat,
        :domain_z                 => opts.domain_z,
        :grid_size                => opts.grid_size,
        :vertical_stretching_mode => opts.vertical_stretching_mode,
        :vertical_grid_file       => opts.vertical_grid_file,
        :vertical_depths          => opts.vertical_depths,
        :bathymetry_source        => Symbol(opts.bathy_source),
        :bathy_source              => opts.bathy_source,
        :coastline_source          => opts.coastline_source,
        :wind_source               => opts.wind_source,
        :wind_time_iso             => opts.wind_time_iso,
        :inshore_depth            => opts.inshore_depth,
        :shelf_slope              => opts.shelf_slope,
        :mask_bay_of_fundy        => opts.mask_bay_of_fundy,
        :atmospheric_source       => opts.atmospheric_source,
        :drag_formulation         => :garratt_1977,
        :bulk_heat_flux           => true,
        :climatology              => false,
        :ocean_boundary_source    => opts.ocean_boundary_source,
        :hydrography_source       => opts.hydrography_source,
        :hydrography_month        => opts.hydrography_month,
        :obc_type                 => opts.obc_type,
        :sponge_layer_width       => opts.sponge_layer_width,
        :sponge_timescale         => opts.sponge_timescale,
        :cd_drag                  => 0.0025,
        :bottom_drag              => 0.0001,
        :bbl_mixing_closure       => :tke,
        :sim_duration_seconds     => opts.sim_duration,
        :sim_dt_seconds           => opts.sim_dt,
        :adaptive_cfl             => opts.adaptive_cfl,
        :target_cfl               => opts.target_cfl,
        :target_wave_cfl          => 0.20,
         :max_dt_seconds           => opts.max_dt,
         # Was hardcoded to 2.0, so the second options struct (built from this dict at the call
         # site below) always carried 2.0 regardless of what the TOML asked for. Use the value
         # that was actually parsed.
         :min_dt_seconds           => opts.min_dt_seconds,
         :coriolis_latitude        => 44.5,
         :divergence_limit         => 20.0,
         :output_path              => opts.output_dir,
         :storage_backend          => opts.storage_backend,
         :output_schedule_seconds  => opts.output_schedule_seconds,
         :enable_checkpoint        => opts.enable_checkpoint,
         :checkpoint_prefix        => opts.checkpoint_prefix,
         :checkpoint_schedule      => opts.checkpoint_schedule,
         :checkpoint_dir           => opts.checkpoint_dir,
         :checkpoint_cleanup       => opts.checkpoint_cleanup,
         :auto_restart             => opts.auto_restart,
         :use_gpu                  => opts.use_gpu,
         :fallback_to_cpu          => opts.fallback_to_cpu
    )
    for (k, v) in kwargs
        d[k] = v
    end

    return HydrodynamicConfig(
        d[:domain_lon], d[:domain_lat], d[:domain_z], d[:grid_size],
        Symbol(d[:vertical_stretching_mode]),         String(d[:vertical_grid_file]),
        Float64[Float64(x) for x in d[:vertical_depths]],
        Symbol(d[:bathymetry_source]),
        String(d[:bathy_source]), String(d[:coastline_source]), String(d[:wind_source]), String(d[:wind_time_iso]), Float64(d[:inshore_depth]),
        Float64(d[:shelf_slope]), Symbol(d[:atmospheric_source]),
        Symbol(d[:drag_formulation]), Bool(d[:bulk_heat_flux]),
        Bool(d[:climatology]),
        Symbol(d[:ocean_boundary_source]), Symbol(d[:obc_type]),
        Float64(d[:sponge_layer_width]), Float64(d[:sponge_timescale]),
        Float64(d[:bottom_drag]), Symbol(d[:bbl_mixing_closure]),
        Float64(d[:sim_duration_seconds]), Float64(d[:sim_dt_seconds]),
        Bool(d[:adaptive_cfl]), Float64(d[:target_cfl]),
        Float64(d[:target_wave_cfl]), Float64(d[:max_dt_seconds]),
        Float64(d[:min_dt_seconds]), Float64(d[:coriolis_latitude]),
         Float64(d[:divergence_limit]), String(d[:output_path]),
         String(d[:storage_backend]), Float64(d[:output_schedule_seconds]),
         Bool(d[:enable_checkpoint]), String(d[:checkpoint_prefix]),
         String(d[:checkpoint_schedule]), String(d[:checkpoint_dir]),
         Bool(d[:checkpoint_cleanup]), Bool(d[:auto_restart]),
         Bool(d[:use_gpu]), Bool(d[:fallback_to_cpu])
    )
end

"""
    to_larval_config(opts::HydrodynamicOptions; kwargs...) -> LarvalDispersalConfig

Extract and convert a `HydrodynamicOptions` instance into a pure `LarvalDispersalConfig`.
"""
function to_larval_config(opts::HydrodynamicOptions; kwargs...)::LarvalDispersalConfig
    d = Dict{Symbol, Any}(
        :n_particles              => opts.n_particles,
        :track_duration_seconds   => opts.track_duration,
        :track_dt_seconds         => opts.track_dt,
        :diffusivity_h            => opts.diffusivity_h,
        :diffusivity_v            => opts.diffusivity_v,
        :min_seabed_depth         => opts.min_seabed_depth,
        :buffer_km                => opts.buffer_km,
        :release_depth_mode       => opts.release_depth_mode,
        :bottom_release_offset    => opts.bottom_release_offset,
        :enable_initial_ascent    => opts.enable_initial_ascent,
        :ascent_speed             => opts.ascent_speed,
        :ascent_target_depth      => opts.ascent_target_depth,
        :enable_dvm               => opts.enable_dvm,
        :zoea1_day_depth          => opts.zoea1_day_depth,
        :zoea1_night_depth        => opts.zoea1_night_depth,
        :zoea2_day_depth          => opts.zoea2_day_depth,
        :zoea2_night_depth        => opts.zoea2_night_depth,
        :megalopa_day_depth       => opts.megalopa_day_depth,
        :megalopa_night_depth     => opts.megalopa_night_depth,
        :zoea2_depth_factor       => opts.zoea2_depth_factor,
        :megalopa_swim_factor     => opts.megalopa_swim_factor,
        :enable_molting           => opts.enable_molting,
        :t_base                   => opts.t_base,
        :dd_zoea1_to_zoea2        => opts.dd_zoea1_to_zoea2,
        :dd_zoea2_to_megalopa     => opts.dd_zoea2_to_megalopa,
        :dd_megalopa_to_settle    => opts.dd_megalopa_to_settle,
        :mortality_base           => opts.mortality_base,
        :mortality_thermal_thresh => opts.mortality_thermal_threshold,
        :mortality_thermal_sens   => opts.mortality_thermal_sensitivity,
        :mortality_cold_thresh    => opts.mortality_cold_threshold,
        :mortality_cold_sens      => opts.mortality_cold_sensitivity,
        :settlement_min_depth     => opts.settlement_min_depth,
        :settlement_max_depth     => opts.settlement_max_depth,
        :settlement_max_temp      => opts.settlement_max_temp,
        :settlement_stochastic    => opts.settlement_stochastic,
        :cv_molt                  => opts.cv_molt,
        :cv_mortality             => opts.cv_mortality,
        :cv_settlement            => opts.cv_settlement,
        :enable_checkpoint        => opts.enable_checkpoint,
        :output_path              => opts.output_dir,
        :storage_backend          => opts.storage_backend,
        :enable_zarr            => opts.enable_zarr,
        :zarr_path              => opts.zarr_path,
        :run_id                   => opts.run_id,
        :seed                     => opts.seed
    )
    for (k, v) in kwargs
        d[k] = v
    end

    return LarvalDispersalConfig(
        Int(d[:n_particles]), Float64(d[:track_duration_seconds]),
        Float64(d[:track_dt_seconds]), Float64(d[:diffusivity_h]),
        Float64(d[:diffusivity_v]), Float64(d[:min_seabed_depth]),
        Float64(d[:buffer_km]), Symbol(d[:release_depth_mode]),
        (Float64(d[:bottom_release_offset][1]), Float64(d[:bottom_release_offset][2])),
        Bool(d[:enable_initial_ascent]), Float64(d[:ascent_speed]),
        Float64(d[:ascent_target_depth]), Bool(d[:enable_dvm]),
        Float64(d[:zoea1_day_depth]), Float64(d[:zoea1_night_depth]),
        Float64(d[:zoea2_day_depth]), Float64(d[:zoea2_night_depth]),
        Float64(d[:megalopa_day_depth]), Float64(d[:megalopa_night_depth]),
        Float64(d[:zoea2_depth_factor]), Float64(d[:megalopa_swim_factor]),
        Bool(d[:enable_molting]), Float64(d[:t_base]),
        Float64(d[:dd_zoea1_to_zoea2]), Float64(d[:dd_zoea2_to_megalopa]),
        Float64(d[:dd_megalopa_to_settle]), Float64(d[:mortality_base]),
        Float64(d[:mortality_thermal_thresh]), Float64(d[:mortality_thermal_sens]),
        Float64(d[:mortality_cold_thresh]), Float64(d[:mortality_cold_sens]),
        Float64(d[:settlement_min_depth]), Float64(d[:settlement_max_depth]),
        Float64(d[:settlement_max_temp]), Bool(d[:settlement_stochastic]),
         Float64(d[:cv_molt]), Float64(d[:cv_mortality]), Float64(d[:cv_settlement]),
         Bool(d[:enable_checkpoint]), String(d[:output_path]),
         String(d[:storage_backend]), Int(d[:seed])
     )
end
