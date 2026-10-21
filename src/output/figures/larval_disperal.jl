"""
    larval_disperal.jl

Larval dispersal figures.

Trajectory maps, vertical-migration profiles, run comparison, settlement
density, empirical movement fields, connectivity and recruitment summaries.
"""

using CairoMakie
import CairoMakie: record

"""
    plot_particle_trajectories(
        trajectories::NamedTuple;
        bathymetry_data::Union{Nothing, NamedTuple} = nothing,
        strata::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        title::AbstractString = "Snow Crab Larval Dispersal Trajectories",
        output_path::Union{Nothing, AbstractString} = "outputs/particle_trajectories.png",
        max_display_particles::Int = 100
    )

Render a 2D spatial map of snow crab larval drift trajectories with optional
background bathymetry and management strata boundaries.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort` with fields
  `(lons, lats, depths, times, ids)`.
- `bathymetry_data::Union{Nothing, NamedTuple}`: Optional `(lon, lat, elevation)`
  for background contours.
- `strata::Union{Nothing, AbstractVector{<:NamedTuple}}`: Optional management strata
  definitions with polygon coordinates.
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Path to save the figure (PNG/PDF).
- `max_display_particles::Int`: Maximum number of individual tracks to render.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- North, E. W., et al. (2009). ICES Cooperative Research Report, No. 295.
"""
function plot_particle_trajectories(
    trajectories::NamedTuple;
    bathymetry_data::Union{Nothing, NamedTuple} = nothing,
    strata::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    title::AbstractString = "Snow Crab Larval Dispersal Trajectories",
    output_path::Union{Nothing, AbstractString} = "outputs/particle_trajectories.png",
    max_display_particles::Int = 100
)
    fig = Figure(size = (950, 720), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    # 1. Plot background bathymetry if provided
    if !isnothing(bathymetry_data)
        co = contourf!(
            ax,
            bathymetry_data.lon,
            bathymetry_data.lat,
            bathymetry_data.elevation,
            colormap = :viridis,
            levels = 15
        )
        Colorbar(fig[1, 2], co, label = "Seafloor Elevation (m)")
    end

    # 2. Draw management strata outlines if provided
    if !isnothing(strata)
        for s in strata
            if hasproperty(s, :polygon) && !isempty(s.polygon)
                poly_lons = [Float64(pt[1]) for pt in s.polygon]
                poly_lats = [Float64(pt[2]) for pt in s.polygon]
                # Close the ring: an unclosed polygon renders as a gap, and a gap in
                # a management boundary reads as an opening that is not there.
                if !isempty(poly_lons) &&
                   (poly_lons[1] != poly_lons[end] || poly_lats[1] != poly_lats[end])
                    push!(poly_lons, poly_lons[1])
                    push!(poly_lats, poly_lats[1])
                end
                lines!(
                    ax,
                    poly_lons,
                    poly_lats,
                    color = (:gold, 0.75),
                    linewidth = 1.8,
                    linestyle = :dash
                )
            end
        end
    end

    # 3. Plot larval tracks with high-visibility outline underlay and vivid color
    n_particles = size(trajectories.lons, 1)
    n_plot = min(n_particles, max_display_particles)

    # High-contrast dark shadow underlay
    for p in 1:n_plot
        lines!(
            ax,
            trajectories.lons[p, :],
            trajectories.lats[p, :],
            color = (:black, 0.75),
            linewidth = 3.8
        )
    end

    # Vivid foreground trajectory lines
    for p in 1:n_plot
        lines!(
            ax,
            trajectories.lons[p, :],
            trajectories.lats[p, :],
            color = :cyan,
            linewidth = 2.2
        )
    end

    # 4. Mark release positions (start) and final drift positions (end) with strokes
    scatter!(
        ax,
        trajectories.lons[1:n_plot, 1],
        trajectories.lats[1:n_plot, 1],
        color = :springgreen,
        strokecolor = :black,
        strokewidth = 1.5,
        marker = :circle,
        markersize = 11,
        label = "Release (t = 0)"
    )
    scatter!(
        ax,
        trajectories.lons[1:n_plot, end],
        trajectories.lats[1:n_plot, end],
        color = :crimson,
        strokecolor = :white,
        strokewidth = 1.5,
        marker = :diamond,
        markersize = 12,
        label = "Final Position"
    )

    axislegend(ax, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_vertical_migration_profiles(
        trajectories::NamedTuple;
        sample_indices = 1:min(5, size(trajectories.depths, 1)),
        title::AbstractString = "Diel Vertical Migration (DVM) Depth Profiles",
        show_surface_target::Bool = true,
        surface_target_depth::Real = -10.0,
        output_path::Union{Nothing, AbstractString} = "outputs/dvm_profiles.png"
    )

Plot larval depth trajectories over time to visualize initial ascent from benthic
release depths and diurnal oscillations between nighttime surface grazing and
daytime sub-surface predator avoidance.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `sample_indices`: Range or vector of particle indices to plot.
- `title::AbstractString`: Figure title.
- `show_surface_target::Bool`: Whether to draw reference line for surface target layer.
- `surface_target_depth::Real`: Surface mixed layer target depth in meters (default -10.0 m).
- `output_path::Union{Nothing, AbstractString}`: Destination path for figure.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- Incze, L. S., et al. (1987). *Marine Biology*, 95(2), 195-200.
- Sainte-Marie, G., & Sainte-Marie, B. (1999). *Can. J. Fish. Aquat. Sci.*, 56(11), 2181-2193.
"""
function plot_vertical_migration_profiles(
    trajectories::NamedTuple;
    sample_indices = 1:min(5, size(trajectories.depths, 1)),
    title::AbstractString = "Diel Vertical Migration (DVM) Depth Profiles",
    show_surface_target::Bool = true,
    surface_target_depth::Real = -10.0,
    output_path::Union{Nothing, AbstractString} = "outputs/dvm_profiles.png"
)
    fig = Figure(size = (850, 480), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Simulation Time (hours)",
        ylabel = "Depth (m)"
    )

    time_hours = trajectories.times ./ 3600.0
    colors = [:royalblue, :darkorange, :forestgreen, :purple, :crimson]

    for (idx, p) in enumerate(sample_indices)
        c = colors[mod1(idx, length(colors))]
        lines!(
            ax,
            time_hours,
            trajectories.depths[p, :],
            color = c,
            linewidth = 2.0,
            label = "Larva $(p)"
        )
        scatter!(
            ax,
            [time_hours[1]],
            [trajectories.depths[p, 1]],
            color = c,
            marker = :circle,
            markersize = 8
        )
    end

    hlines!(ax, [0.0], color = :gray50, linestyle = :dash)
    if show_surface_target
        hlines!(
            ax,
            [Float64(surface_target_depth)],
            color = :teal,
            linestyle = :dot,
            linewidth = 1.5,
            label = "Surface Target ($(round(surface_target_depth, digits=1))m)"
        )
    end
    axislegend(ax, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    compare_run_dispersal(
        run_dict::Union{AbstractDict, NamedTuple};
        title::AbstractString = "Larval Dispersal Across Runs",
        output_path::Union{Nothing, AbstractString} = "outputs/run_comparison.png"
    )

Compare spatial trajectories and final dispersion centroids across several
archived runs.

# Mathematical Formulations
- **Centroid Center of Mass**:
  ```math
  \\bar{\\\\\\\\lambda} = \\\\\\\\\\\\\\frac{1}{N_p} \\sum_{p=1}^{N_p} \\\\\\\\lambda_p(t_{\\text{end}}), \\quad
  \\bar{\\phi} = \\\\\\\\\\\\\\frac{1}{N_p} \\sum_{p=1}^{N_p} \\phi_p(t_{\\text{end}})
  ```

# Inputs
- `run_dict::Union{AbstractDict, NamedTuple}`: Mapping of run identifiers or keys
  to trajectory NamedTuples (or data containers holding a `.trajectories` field).
- `title::AbstractString`: Comparison figure title.
- `output_path::Union{Nothing, AbstractString}`: Destination filepath to save output.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- Brickman, D., et al. (2018). *Progress in Oceanography*, 164, 49-64.
"""
function compare_run_dispersal(
    run_dict::Union{AbstractDict, NamedTuple};
    title::AbstractString = "Larval Dispersal Across Runs",
    output_path::Union{Nothing, AbstractString} = "outputs/run_comparison.png"
)
    n_runs = length(run_dict)
    n_cols = min(n_runs, 3)
    n_rows = ceil(Int, n_runs / 3)
    fig = Figure(size = (450 * n_cols, 420 * n_rows), fontsize = 12)

    run_keys = run_dict isa NamedTuple ? keys(run_dict) : collect(keys(run_dict))
    palette = (:royalblue, :seagreen, :darkorange, :orangered, :crimson, :firebrick)

    for (r_idx, r_key) in enumerate(run_keys)
        row = div(r_idx - 1, 3) + 1
        col = mod(r_idx - 1, 3) + 1

        raw_val = run_dict isa NamedTuple ? getproperty(run_dict, r_key) : run_dict[r_key]
        traj = hasproperty(raw_val, :trajectories) ? raw_val.trajectories : raw_val
        c = palette[mod1(r_idx, length(palette))]

        ax = Axis(
            fig[row, col],
            title = "$(string(r_key))",
            xlabel = "Longitude (°E)",
            ylabel = "Latitude (°N)"
        )

        n_p = size(traj.lons, 1)
        for p in 1:min(n_p, 50)
            lines!(
                ax, traj.lons[p, :], traj.lats[p, :],
                color = (:black, 0.5), linewidth = 2.8
            )
            lines!(
                ax, traj.lons[p, :], traj.lats[p, :],
                color = (c, 0.85), linewidth = 1.8
            )
        end

        mean_lon_end = sum(traj.lons[:, end]) / n_p
        mean_lat_end = sum(traj.lats[:, end]) / n_p

        scatter!(
            ax,
            [mean_lon_end],
            [mean_lat_end],
            color = :gold,
            strokecolor = :black,
            strokewidth = 1.5,
            marker = :star5,
            markersize = 16,
            label = "Centroid"
        )
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_larval_dispersal_density(
        trajectories::NamedTuple;
        lon_bins = range(-68.0, -57.0, length = 100),
        lat_bins = range(41.0, 48.0, length = 100),
        title::AbstractString = "Larval Settlement Density Kernel",
        output_path::Union{Nothing, AbstractString} = "outputs/settlement_density.png"
    )

Compute and plot a 2D histogram heatmap of final particle positions representing
larval settlement and nursery retention hotspots.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `lon_bins`: Grid longitude bin edges.
- `lat_bins`: Grid latitude bin edges.
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Path to save figure.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_larval_dispersal_density(
    trajectories::NamedTuple;
    lon_bins = range(-68.0, -57.0, length = 100),
    lat_bins = range(41.0, 48.0, length = 100),
    title::AbstractString = "Larval Settlement Density Kernel",
    output_path::Union{Nothing, AbstractString} = "outputs/settlement_density.png"
)
    fig = Figure(size = (850, 650), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    end_lons = trajectories.lons[:, end]
    end_lats = trajectories.lats[:, end]
    n_p = length(end_lons)

    nx = length(lon_bins) - 1
    ny = length(lat_bins) - 1
    density = zeros(Float64, nx, ny)

    for p in 1:n_p
        x = end_lons[p]
        y = end_lats[p]
        if lon_bins[1] <= x <= lon_bins[end] && lat_bins[1] <= y <= lat_bins[end]
            i = clamp(searchsortedlast(lon_bins, x), 1, nx)
            j = clamp(searchsortedlast(lat_bins, y), 1, ny)
            density[i, j] += 1.0
        end
    end

    density ./= max(1.0, sum(density))
    density .*= 100.0

    lon_centers = [(lon_bins[i] + lon_bins[i + 1]) / 2.0 for i in 1:nx]
    lat_centers = [(lat_bins[j] + lat_bins[j + 1]) / 2.0 for j in 1:ny]

    hm = heatmap!(ax, lon_centers, lat_centers, density, colormap = :inferno)
    Colorbar(fig[1, 2], hm, label = "Settlement Density (%)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_empirical_movement_field(
        emp_mov::NamedTuple;
        title::AbstractString = "Empirical Advection & Turbulent Diffusivity Field",
        output_path::Union{Nothing, AbstractString} = "outputs/empirical_movement.png"
    )

Plot empirical velocity vector arrows over the empirical turbulent diffusivity field.

# Inputs
- `emp_mov::NamedTuple`: Output from `estimate_empirical_movement`.
- `title::AbstractString`: Plot title.
- `output_path::Union{Nothing, AbstractString}`: File path to save figure.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_empirical_movement_field(
    emp_mov::NamedTuple;
    title::AbstractString = "Empirical Advection & Turbulent Diffusivity Field",
    output_path::Union{Nothing, AbstractString} = "outputs/empirical_movement.png"
)
    fig = Figure(size = (850, 650), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )

    lon_c = emp_mov.lon_centers
    lat_c = emp_mov.lat_centers
    diff = copy(emp_mov.diffusivity)
    diff[isnan.(diff)] .= 0.0

    hm = heatmap!(ax, lon_c, lat_c, diff, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Empirical Diffusivity D (m² s⁻¹)")

    pts = Point2f[]
    dirs = Vec2f[]
    u_m = emp_mov.u_mean
    v_m = emp_mov.v_mean

    for i in 1:length(lon_c), j in 1:length(lat_c)
        if !isnan(u_m[i, j]) && !isnan(v_m[i, j]) && emp_mov.sample_count[i, j] > 2
            push!(pts, Point2f(lon_c[i], lat_c[j]))
            push!(dirs, Vec2f(u_m[i, j] * 0.5, v_m[i, j] * 0.5))
        end
    end

    if !isempty(pts)
        if isdefined(Makie, :arrows2d!)
            arrows2d!(
                ax,
                pts,
                dirs,
                tipwidth = 6,
                tiplength = 10,
                lengthscale = 1.0,
                tipcolor = :white,
                shaftcolor = :white
            )
        else
            arrows!(
                ax,
                pts,
                dirs,
                tipwidth = 6,
                tiplength = 10,
                lengthscale = 1.0,
                tipcolor = :white,
                shaftcolor = :white
            )
        end
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_connectivity_matrix(
        conn::NamedTuple;
        title::AbstractString = "Macro-Regional Population Connectivity Matrix",
        output_path::Union{Nothing, AbstractString} = "outputs/connectivity_matrix.png"
    )

Render an annotated transition probability heatmap showing larval transfer probabilities
and self-retention percentages between management zones.

# Inputs
- `conn::NamedTuple`: Output from `compute_empirical_connectivity`.
- `title::AbstractString`: Title string.
- `output_path::Union{Nothing, AbstractString}`: Path for output image.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_connectivity_matrix(
    conn::NamedTuple;
    title::AbstractString = "Macro-Regional Population Connectivity Matrix",
    output_path::Union{Nothing, AbstractString} = "outputs/connectivity_matrix.png"
)
    mat = conn.matrix .* 100.0
    strata = conn.strata_names
    n_s = length(strata)

    fig = Figure(size = (800, 650), fontsize = 12)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Destination Stratum",
        ylabel = "Source Stratum",
        xticks = (1:n_s, strata),
        yticks = (1:n_s, strata),
        xticklabelrotation = 0.35
    )

    hm = heatmap!(ax, 1:n_s, 1:n_s, mat', colormap = :blues)
    Colorbar(fig[1, 2], hm, label = "Transition Probability (%)")

    for i in 1:n_s, j in 1:n_s
        val = mat[i, j]
        txt_col = val > 50.0 ? :white : :black
        text!(
            ax,
            Point2f(j, i),
            text = string(round(val, digits = 1), "%"),
            align = (:center, :center),
            color = txt_col,
            fontsize = 12
        )
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_thermal_exposure_map(
        therm_met::NamedTuple;
        title::AbstractString = "Larval Thermal Exposure & Degree-Days",
        output_path::Union{Nothing, AbstractString} = "outputs/thermal_exposure.png"
    )

Plot a two-panel map showing accumulated thermal degree-days and mean exposure temperature.
"""
function plot_thermal_exposure_map(
    therm_met::NamedTuple;
    title::AbstractString = "Larval Thermal Exposure & Degree-Days",
    output_path::Union{Nothing, AbstractString} = "outputs/thermal_exposure.png"
)
    fig = Figure(size = (1100, 500), fontsize = 12)

    lon_c = therm_met.lon_centers
    lat_c = therm_met.lat_centers

    # Panel 1: Accumulated Degree-Days
    ax1 = Axis(
        fig[1, 1],
        title = "Mean Cumulative Degree-Days (°C·days)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    dd = copy(therm_met.mean_degree_days)
    dd[isnan.(dd)] .= 0.0
    hm1 = heatmap!(ax1, lon_c, lat_c, dd, colormap = :inferno)
    Colorbar(fig[1, 2], hm1, label = "Degree-Days (°C·d)")

    # Panel 2: Mean Ambient Exposure Temperature
    ax2 = Axis(
        fig[1, 3],
        title = "Mean Exposure Temperature (°C)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    temp = copy(therm_met.mean_exposure_temperature)
    temp[isnan.(temp)] .= 0.0
    hm2 = heatmap!(ax2, lon_c, lat_c, temp, colormap = :thermal)
    Colorbar(fig[1, 4], hm2, label = "Temperature (°C)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_recruitment_summary(
        rec_met::NamedTuple;
        title::AbstractString = "Larval Recruitment & Retention Summary",
        output_path::Union{Nothing, AbstractString} = "outputs/recruitment_summary.png"
    )

Plot a multi-panel spatial distribution of spawning release, settlement density, and recruitment success rate.
"""
function plot_recruitment_summary(
    rec_met::NamedTuple;
    title::AbstractString = "Larval Recruitment & Retention Summary",
    output_path::Union{Nothing, AbstractString} = "outputs/recruitment_summary.png"
)
    fig = Figure(size = (1100, 500), fontsize = 12)
    lon_c = rec_met.lon_centers
    lat_c = rec_met.lat_centers

    # Panel 1: Settlement Density
    ax1 = Axis(
        fig[1, 1],
        title = "Settlement Density (count)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    hm1 = heatmap!(ax1, lon_c, lat_c, rec_met.settlement_density, colormap = :plasma)
    Colorbar(fig[1, 2], hm1, label = "Settled Larvae")

    # Panel 2: Recruitment Success Rate
    ax2 = Axis(
        fig[1, 3],
        title = "Benthic Recruitment Success Rate (%)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)"
    )
    succ_pct = rec_met.success_rate .* 100.0
    hm2 = heatmap!(ax2, lon_c, lat_c, succ_pct, colormap = :viridis)
    Colorbar(fig[1, 4], hm2, label = "Success Rate (%)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    _timeseries_value(frame::NamedTuple, i::Int, j::Int, k::Int, variable::Symbol) -> Float64

Read one field of an extracted hydrodynamic frame at grid point `(i, j, k)`, normalising the
scalars that the plots conventionally rescale (Brunt-Väisälä, diffusivity and speed are drawn in
their scientific units of 10⁻⁴ and cm s⁻¹).
"""
@inline function _timeseries_value(frame::NamedTuple, i::Int, j::Int, k::Int, variable::Symbol)
    v = if variable === :salinity || variable === :S
        frame.salinity[i, j, k]
    elseif variable === :stratification || variable === :N2
        frame.stratification[i, j, k] * 1.0e4
    elseif variable === :diffusion || variable === :kappa
        frame.diffusion[i, j, k] * 1.0e4
    elseif variable === :viscosity || variable === :nu
        frame.viscosity[i, j, k] * 1.0e4
    elseif variable === :speed
        frame.speed[i, j, k] * 100.0
    elseif variable === :density
        frame.density[i, j, k]
    elseif variable === :vorticity
        frame.vorticity[i, j, k]
    elseif variable === :w
        frame.w[i, j, k]
    else
        frame.temperature[i, j, k]
    end
    return Float64(v)
end

"""
    plot_hydrodynamic_hovmoller(
        hydrodynamics::Any;
        station = nothing, lon = -63.5, lat = 44.0,
        variable = :temperature, max_frames = 240,
        title = nothing, output_path = "outputs/hydrodynamic_hovmoller.png"
    ) -> Figure

Time–depth (Hovmöller) diagram of a hydrodynamic variable at a fixed station, built from the real
snapshots.

This is the canonical Eulerian diagnostic for a shelf model and the single most useful missing
figure: it shows at once whether the simulation resolves the structures that matter here — the cold
intermediate layer as a sloping isotherm band, the M2/S2 tide as horizontal banding, and any
periodic or drifting behaviour. A time series alone cannot show which of those is present, because
it collapses the vertical structure.

# Inputs
- `hydrodynamics`: Model output, or a path the extractor understands.
- `station`: `(lon, lat)` tuple, or `nothing` to use `lon`/`lat`.
- `variable`: `:temperature`, `:salinity`, `:speed`, `:w`, or `:density`.
- `max_frames`: Cap on snapshots walked; they are evenly subsampled beyond this.
- `output_path`: PNG path, or `nothing`.
"""
function plot_hydrodynamic_hovmoller(
    hydrodynamics::Any;
    station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    lon::Real = -63.5,
    lat::Real = 44.0,
    variable::Symbol = :temperature,
    max_frames::Int = 240,
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_hovmoller.png"
)
    first_frame = extract_hydrodynamic_dataset(hydrodynamics)
    target_lon = !isnothing(station) ? Float64(station[1]) : Float64(lon)
    target_lat = !isnothing(station) ? Float64(station[2]) : Float64(lat)
    i_st = argmin(abs.(first_frame.lons .- target_lon))
    j_st = argmin(abs.(first_frame.lats .- target_lat))

    n_snap = isdefined(first_frame, :n_snapshots) ? first_frame.n_snapshots : length(first_frame.times)
    idx = n_snap <= max_frames ? collect(1:n_snap) :
          unique(round.(Int, range(1, n_snap, length = max_frames)))
    length(idx) < 2 && (idx = [1, max(1, n_snap)])

    zs = collect(Float64, first_frame.depths)
    nz = length(zs)
    M = fill(NaN, nz, length(idx))
    hours = Float64[]
    for (c, k) in enumerate(idx)
        f = k == 1 ? first_frame :
            extract_hydrodynamic_dataset(hydrodynamics; time_index = k)
        push!(hours, first(f.times) / 3600.0)
        for r in 1:nz
            M[r, c] = _timeseries_value(f, i_st, j_st, r, variable)
        end
    end

    var_title = uppercase(string(variable)[1:1]) * string(variable)[2:end]
    fig = Figure(size = (980, 520), fontsize = 12)
    ax = Axis(fig[1, 1],
              title = something(title, "Time–Depth ($(var_title)) at " *
                                      "($(round(first_frame.lons[i_st], digits=2))°E, " *
                                      "$(round(first_frame.lats[j_st], digits=2))°N)"),
              xlabel = "Simulation time (hours)", ylabel = "Depth (m)")
    hm = heatmap!(ax, hours, zs, M; colormap = :balance,
                   colorrange = _field_colorrange(M))
    Colorbar(fig[1, 2], hm, label = var_title)
    # Isobaths: shading the T–S–z space without them hides the structure this plot exists to show.
    isotherms = variable === :temperature
    isotherms && vlines!(ax, [12.42, 24.0], color = (:white, 0.35), linestyle = :dot,
                         label = "M2 / S2")
    isotherms && axislegend(ax, position = :rb)
    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end

"""
    plot_bottom_temperature_map(hydrodynamics; depth = nothing, output_path = ...)

Map of near-bed temperature across the domain.

Settlement suitability in the larval model is driven by near-bed temperature, and the recruitment
figures report *where* larvae settled but never the temperature field that decided it. Mapping it
next to the settlement-density figure makes it possible to see whether a recruitment pattern
follows bottom temperature or contradicts it.

Land is masked to NaN so the map shows the marine domain only.
"""
function plot_bottom_temperature_map(
    hydrodynamics::Any;
    depth::Union{Nothing, Real} = nothing,
    domain_lon::Tuple{<:Real, <:Real} = (-71.0, -53.0),
    domain_lat::Tuple{<:Real, <:Real} = (40.0, 48.5),
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/bottom_temperature_map.png"
)
    f = extract_hydrodynamic_dataset(hydrodynamics; depth = depth,
                                    domain_lon = domain_lon, domain_lat = domain_lat)
    zs = collect(Float64, f.depths)
    k = isnothing(depth) ? argmax(zs) : resolve_depth_index(zs, depth, nothing)  # deepest available
    Tbot = permutedims(f.temperature[:, :, k])
    fig = Figure(size = (760, 520), fontsize = 12)
    ax = Axis(fig[1, 1],
              title = something(title, "Near-Bed Temperature (deepest resolved level, $(zs[k]) m)"),
              xlabel = "Longitude (°E)", ylabel = "Latitude (°N)")
    hm = heatmap!(ax, f.lons, f.lats, Tbot, colormap = :thermal;
                   colorrange = _field_colorrange(Tbot))
    Colorbar(fig[1, 2], hm, label = "Temperature (°C)")
    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end

"""
    plot_temperature_salinity_diagram(hydrodynamics; output_path = ...)

Temperature–salinity (T–S) scatter of the domain, coloured by depth.

The standard water-mass diagnostic. For this domain it answers whether the simulated water column
contains the three masses the larvae and the shelf exchange depend on — cold fresh shelf water
(Gulf of St. Lawrence), the cold intermediate layer, and warm salty slope water — and whether the
model produces a continuous mixing line between them or an implausibly disconnected set.
"""
function plot_temperature_salinity_diagram(
    hydrodynamics::Any;
    domain_lon::Tuple{<:Real, <:Real} = (-71.0, -53.0),
    domain_lat::Tuple{<:Real, <:Real} = (40.0, 48.5),
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/ts_diagram.png"
)
    f = extract_hydrodynamic_dataset(hydrodynamics; domain_lon = domain_lon, domain_lat = domain_lat)
    zs = collect(Float64, f.depths)
    # Scatter every (lon, lat, depth) cell, so vectorise the 3-D fields rather than reducing.
    Sv = vec(f.salinity)
    Tv = vec(f.temperature)
    Zv = repeat(zs, outer = div(length(Sv), length(zs)))

    fig = Figure(size = (760, 560), fontsize = 12)
    ax = Axis(fig[1, 1], title = something(title, "Temperature–Salinity Diagram"),
              xlabel = "Salinity (psu)", ylabel = "Temperature (°C)")
    m = isfinite.(Tv) .& isfinite.(Sv)
    sc = scatter!(ax, Sv[m], Tv[m]; color = Zv[m], colormap = :viridis)
    Colorbar(fig[1, 2], sc, label = "Depth (m)")
    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end

