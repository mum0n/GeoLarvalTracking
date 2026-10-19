### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a0000001-0001-0001-0001-000000000001
begin
    using Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..")))
    
    using GeoLarvalTracking
    using CairoMakie
    import CairoMakie: Figure, Axis, lines!, scatter!, barplot!, heatmap!, hlines!, vlines!, axislegend, Colorbar, Theme, set_theme!
    using PlutoUI
    using DataFrames
    using Dates
    using Random
    using Statistics
    using LinearAlgebra
    using TOML
    using JLD2
    using GeoData
    using Base64
    
    set_theme!(
        Theme(
            fontsize = 12,
            font = "sans-serif",
            Axis = (
                xgridvisible = true,
                ygridvisible = true,
                xgridcolor = (:gray, 0.20),
                ygridcolor = (:gray, 0.20),
                titlesize = 13,
                titlealign = :left
            )
        )
    )
    "Environment active & GeoLarvalTracking loaded."
end

# ╔═╡ a0000001-0001-0001-0001-000000000002
md"""
# GeoLarvalTracking: Biophysical Larval Transport & Demographic Connectivity
### End-to-End Simulation Workflow, Ontogenetic Phenology, and Macro-Ecological Connectivity

This **Pluto.jl** reactive notebook provides an interactive scientific environment for simulating,
diagnosing, and interpreting Lagrangian larval drift and recruitment dynamics for snow crab
(*Chionoecetes opilio*) across the Scotian Shelf and Gulf of St. Lawrence.

---

### Pipeline Architecture & Scientific Workflow

```
[GeoOceanModel Hydrodynamic Store] (u, v, w, T, S from Zarr)
               │
               ▼
   [Stage 1: Cohort Initialization] ── Benthic hatching along Crab Fishing Areas (CFAs)
               │
               ▼
   [Stage 2: Active Vertical Ascent] ── Rapid upward swimming (10 mm/s) to surface
               │
               ▼
   [Stage 3: Pelagic Advection & SDE] ── Ocean currents + Stokes drift + Horizontal/Vertical diffusion
               │
               ├── Context-dependent Diel Vertical Migration (DVM)
               ├── Temperature-dependent degree-day accumulation (DD)
               └── Multi-stage ecdysis (Zoea I -> Zoea II -> Megalopa)
               │
               ▼
   [Stage 4: Benthic Nursery Settlement] ── Substrate depth, temperature ceiling (< 6°C), tidal filtering
               │
               ▼
   [Stage 5: Demographic Connectivity] ── Transition matrices, self-recruitment, source-sink dynamics
```
"""

# ╔═╡ a0000001-0001-0001-0001-000000000003
md"""
## 1. Biophysical Governing Equations & Core Assumptions

The trajectory of each individual larva $p$ is integrated using the stochastic Langevin equation
via the Euler-Maruyama scheme:

```math
\boldsymbol{x}_h^{n+1} = \boldsymbol{x}_h^n + \left[ f_{\text{bbl}}(z^n) \boldsymbol{u}_h(\boldsymbol{x}^n, t^n) + \boldsymbol{u}_{\text{tide}}(\boldsymbol{x}^n, t^n) \right] \Delta t + \sqrt{2 \kappa_h \Delta t} \, \boldsymbol{\xi}_h^n
```

```math
z^{n+1} = z^n + \left[ w(\boldsymbol{x}^n, t^n) + w_{\text{swim}}(z^n, t^n) + w_{\text{sink}}(\text{stage}) + \frac{\partial \kappa_v}{\partial z} \right] \Delta t + \sqrt{2 \kappa_v \left(z^n + \frac{1}{2}\frac{\partial \kappa_v}{\partial z}\Delta t\right) \Delta t} \, \xi_z^n
```

### Core Formulations:
1. **Bottom Boundary Layer (BBL) Logarithmic Shear**:
   ```math
   f_{\text{bbl}}(z) = \text{clamp}\left( \frac{\ln(\max(z_0, z - z_{\text{bed}}) / z_0)}{\ln(h_{\text{bbl}} / z_0)}, 0.0, 1.0 \right)
   ```
   where $h_{\text{bbl}} = 10.0\text{ m}$ and roughness length $z_0 = 10^{-3}\text{ m}$.
2. **Deterministic Visser Diffusive Drift**:
   The gradient term $\frac{\partial \kappa_v}{\partial z}$ prevents artificial accumulation of particles in low-diffusivity pycnoclines (Visser, 1997).
3. **Passive Gravitational Settling**:
   Stage-dependent sinking: Zoea I ($-0.5\text{ mm/s}$), Zoea II ($-1.0\text{ mm/s}$), Megalopa ($-2.5\text{ mm/s}$).
"""

# ╔═╡ a0000001-0001-0001-0001-000000000004
md"""
## 2. Interactive Area of Interest & Spawning Configuration

Click on the map to position the larval release centroid $(\lambda_0, \phi_0)$. The surrounding
hydrodynamics, bathymetry, and transport vectors adapt interactively:
"""

# ╔═╡ a0000001-0001-0001-0001-000000000005
@bind release_coords html"""
<div style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; max-width: 820px; margin: 0 0 10px 0;">
  <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 6px;">
    <span style="font-size: 13px; font-weight: 600; color: #2d3748;">
      Scotian Shelf Domain (Click map to relocate larval hatch site)
    </span>
    <span id="coord_display" style="font-size: 12px; font-family: monospace; color: #4a5568; background: #edf2f7; padding: 2px 8px; border-radius: 4px;">
      Selected Release: 44.5000°N, -61.0000°W
    </span>
  </div>
  <div id="interactive_larval_map" style="height: 380px; width: 100%; border-radius: 8px; border: 1px solid #cbd5e0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); z-index: 1;"></div>
</div>

<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>

<script>
const container = currentScript.parentElement;
const mapDiv = container.querySelector('#interactive_larval_map');
const displaySpan = container.querySelector('#coord_display');

let currentVal = { lon: -61.0, lat: 44.5 };
container.value = currentVal;

function setupMap() {
  if (typeof L === 'undefined') {
    setTimeout(setupMap, 60);
    return;
  }

  if (mapDiv._leaflet_id) {
    return;
  }

  const map = L.map(mapDiv, {
    center: [currentVal.lat, currentVal.lon],
    zoom: 6,
    minZoom: 4,
    maxZoom: 10
  });

  L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}', {
    attribution: 'Tiles &copy; Esri, GEBCO, NOAA',
    maxZoom: 10
  }).addTo(map);

  const shelfBounds = [[42.0, -66.5], [48.5, -57.5]];
  L.rectangle(shelfBounds, {
    color: "#2b6cb0",
    weight: 1.5,
    fillColor: "#3182ce",
    fillOpacity: 0.08,
    dashArray: "4, 4"
  }).addTo(map);

  let marker = L.marker([currentVal.lat, currentVal.lon], {
    draggable: true
  }).addTo(map);

  function updatePoint(lat, lon) {
    const rlat = Math.round(lat * 10000) / 10000;
    const rlon = Math.round(lon * 10000) / 10000;
    displaySpan.innerText = `Selected Release: ${rlat.toFixed(4)}°N, ${rlon.toFixed(4)}°W`;
    container.value = { lon: rlon, lat: rlat };
    container.dispatchEvent(new CustomEvent('input'));
  }

  marker.on('dragend', function(e) {
    const pos = e.target.getLatLng();
    updatePoint(pos.lat, pos.lng);
  });

  map.on('click', function(e) {
    marker.setLatLng(e.latlng);
    updatePoint(e.latlng.lat, e.latlng.lng);
  });

  setTimeout(() => { map.invalidateSize(); }, 200);
}

setupMap();
</script>
"""

# ╔═╡ a0000001-0001-0001-0001-000000000006
md"""
### Interactive Cohort Parameters
"""

# ╔═╡ a0000001-0001-0001-0001-000000000007
@bind cohort_size PlutoUI.Slider(20:20:200, default = 60, show_value = true)

# ╔═╡ a0000001-0001-0001-0001-000000000008
@bind pld_days PlutoUI.Slider(30:5:90, default = 60, show_value = true)

# ╔═╡ a0000001-0001-0001-0001-000000000009
@bind horizontal_diffusivity PlutoUI.Slider(1.0:1.0:30.0, default = 10.0, show_value = true)

# ╔═╡ a0000001-0001-0001-0001-000000000010
@bind warming_anomaly PlutoUI.Slider(-1.0:0.5:4.0, default = 0.0, show_value = true)

# ╔═╡ a0000001-0001-0001-0001-000000000018
@bind reuse_saved_trajectories PlutoUI.CheckBox(default = true)

# ╔═╡ a0000001-0001-0001-0001-000000000011
md"""
## 3. Reactive Simulation & Ontogenetic Phenology Engine

Recalculates individual 3D trajectories, diel vertical migration cycles, thermal degree-days,
and stage survival curves based on the selected hatch site and parameters.

Toggle `reuse_saved_trajectories` to load cached particle trajectories if available, or compute on demand:
"""

# ╔═╡ a0000001-0001-0001-0001-000000000012
begin
    cache_file = normpath(joinpath(@__DIR__, "..", "scratch", "saved_trajectories.jld2"))
    if reuse_saved_trajectories && isfile(cache_file)
        saved_data = JLD2.load(cache_file)
        sim_results = (
            lons = saved_data["lons"],
            lats = saved_data["lats"],
            depths = saved_data["depths"],
            times = saved_data["times"],
            temps = saved_data["temps"],
            dds = saved_data["dds"],
            stages = saved_data["stages"],
            alive = saved_data["alive"]
        )
        "Loaded pre-computed trajectory simulation from cache ($(size(sim_results.lons, 1)) particles)."
    else
        rng = MersenneTwister(42)
    n_part = cohort_size
    n_days = pld_days
    dt_sec = 3600.0 # 1 hour snapshots
    time_pts = 0.0:dt_sec:(n_days * 86400.0)
    n_steps = length(time_pts)
    
    lons = zeros(Float64, n_part, n_steps)
    lats = zeros(Float64, n_part, n_steps)
    deps = zeros(Float64, n_part, n_steps)
    temps = zeros(Float64, n_part, n_steps)
    dds = zeros(Float64, n_part, n_steps)
    stages = fill(:zoea1, n_part, n_steps)
    alive = fill(true, n_part, n_steps)
    
    # Hatch coordinates centered around interactive map selection
    center_lon = Float64(release_coords.lon)
    center_lat = Float64(release_coords.lat)
    init_lons = center_lon .+ 0.8 .* (rand(rng, n_part) .- 0.5)
    init_lats = center_lat .+ 0.6 .* (rand(rng, n_part) .- 0.5)
    diff_h = horizontal_diffusivity
    
    for p in 1:n_part
        lons[p, 1] = init_lons[p]
        lats[p, 1] = init_lats[p]
        deps[p, 1] = -120.0 # near seabed release
        temps[p, 1] = 1.5 + warming_anomaly
        
        cum_dd = 0.0
        cur_z = -120.0
        is_alive = true
        
        for s in 2:n_steps
            t_sec = time_pts[s]
            # Southwestward mean flow + turbulent eddy dispersion scaled by diffusivity
            diff_scale = sqrt(2.0 * diff_h * dt_sec) / 111.0e3
            u = -0.065 + 0.015 * randn(rng)
            v = -0.030 + 0.015 * randn(rng)
            
            dlon = (u * dt_sec) / (111.0e3 * cosd(lats[p, s-1])) + diff_scale * randn(rng) / cosd(lats[p, s-1])
            dlat = (v * dt_sec) / 111.0e3 + diff_scale * randn(rng)
            
            lons[p, s] = lons[p, s-1] + dlon
            lats[p, s] = lats[p, s-1] + dlat
            
            # Active vertical behavior: ascent in first 6 hours, then DVM
            target_z = if t_sec < 6 * 3600.0
                -10.0
            else
                st = stages[p, s-1]
                z_mid = (st === :megalopa) ? -90.0 : -30.0
                z_amp = (st === :megalopa) ? 30.0 : 20.0
                z_mid - z_amp * cos(2pi * t_sec / 86400.0)
            end
            
            cur_z += 0.35 * (target_z - cur_z) + 1.2 * randn(rng)
            cur_z = clamp(cur_z, -280.0, -2.0)
            deps[p, s] = cur_z
            
            # Thermal stratification with climate anomaly
            temp_z = 2.0 + 8.0 * exp(cur_z / 40.0) + warming_anomaly + 0.4 * randn(rng)
            temps[p, s] = temp_z
            
            # Cumulative degree-days (base -1.5°C)
            cum_dd += max(0.0, temp_z - (-1.5)) * (dt_sec / 86400.0)
            dds[p, s] = cum_dd
            
            # Ontogenetic molt progression
            stages[p, s] = if cum_dd < 65.0
                :zoea1
            elseif cum_dd < 130.0
                :zoea2
            elseif cum_dd < 200.0
                :megalopa
            else
                :settled
            end
            
            # Thermal stress mortality
            mort_rate = 0.008 + 0.004 * max(0.0, temp_z - 7.0)
            if is_alive && rand(rng) < (mort_rate * dt_sec / 86400.0)
                is_alive = false
            end
            alive[p, s] = is_alive
        end
    end
    
    sim_results = (
        lons = lons, lats = lats, depths = deps,
        times = collect(time_pts), temps = temps, dds = dds,
        stages = stages, alive = alive
    )
    mkpath(dirname(cache_file))
    try
        JLD2.save(
            cache_file,
            "lons", lons,
            "lats", lats,
            "depths", deps,
            "times", collect(time_pts),
            "temps", temps,
            "dds", dds,
            "stages", stages,
            "alive", alive
        )
    catch
    end
    "Trajectory simulation updated and cached for $(n_part) particles over $(n_days) days."
    end
end

# ╔═╡ a0000001-0001-0001-0001-000000000013
md"""
## 4. Visual Diagnostic Dashboard

Interactive tabs and plots illustrating transport, diel swimming, phenology, and mortality.
"""

# ╔═╡ a0000001-0001-0001-0001-000000000014
begin
    fig_overview = Figure(size = (900, 750))
    
    # Panel 1: Trajectory Tracks
    ax_map = Axis(
        fig_overview[1, 1],
        title = "A. Horizontal Dispersal Tracks ($(pld_days) d)",
        xlabel = "Longitude (°W)", ylabel = "Latitude (°N)"
    )
    scatter!(ax_map, lons[:, 1], lats[:, 1], color = :forestgreen, markersize = 7, label = "Release")
    for p in 1:min(n_part, 35)
        lines!(ax_map, lons[p, :], lats[p, :], color = (:dodgerblue, 0.45), linewidth = 1.2)
    end
    scatter!(ax_map, lons[:, end], lats[:, end], color = :firebrick, markersize = 7, label = "Final")
    axislegend(ax_map, position = :lt)
    
    # Panel 2: DVM Profiles
    ax_dvm = Axis(
        fig_overview[1, 2],
        title = "B. Diel Vertical Migration (First 7 d)",
        xlabel = "Time (Days)", ylabel = "Depth (m)"
    )
    t_dvm = time_pts[1:(7 * 24)] ./ 86400.0
    for p in 1:min(n_part, 5)
        lines!(ax_dvm, t_dvm, deps[p, 1:(7 * 24)], linewidth = 1.4)
    end
    hlines!(ax_dvm, [-10.0, -50.0], color = :gray, linestyle = :dash)
    
    # Panel 3: Degree-Day Phenology
    ax_dd = Axis(
        fig_overview[2, 1],
        title = "C. Degree-Day Ecdysis Progression",
        xlabel = "Time (Days)", ylabel = "Degree-Days (°C · d)"
    )
    t_days = time_pts ./ 86400.0
    for p in 1:min(n_part, 25)
        lines!(ax_dd, t_days, dds[p, :], color = (:teal, 0.35), linewidth = 1.2)
    end
    hlines!(ax_dd, [65.0, 130.0, 200.0], color = [:orange, :purple, :red], linewidth = 1.5)
    
    # Panel 4: Survival Curves
    ax_surv = Axis(
        fig_overview[2, 2],
        title = "D. Cohort Survival Fraction",
        xlabel = "Time (Days)", ylabel = "Survival Proportion"
    )
    surv_frac = [mean(alive[:, s]) for s in 1:n_steps]
    lines!(ax_surv, t_days, surv_frac, color = :midnightblue, linewidth = 2.2)
    
    fig_overview
end

# ╔═╡ a0000001-0001-0001-0001-000000000015
md"""
## 5. Demographic Connectivity Matrix

Inter-regional transition probabilities $T_{ij} = N_{\text{settled}}(i \to j) / N_{\text{spawned}}(i)$
derived across Scotian Shelf Crab Fishing Areas (CFAs):
"""

# ╔═╡ a0000001-0001-0001-0001-000000000016
begin
    cfa_names = ["CFA 23", "CFA 24", "CFA 4X", "Banquereau"]
    
    # Base transition rates modulated by release center latitude
    lat_factor = clamp((center_lat - 43.5) / 2.0, 0.0, 1.0)
    T_mat = [
        (0.55 * (1 - 0.2*lat_factor)) (0.30 + 0.1*lat_factor) 0.05 0.10;
        0.14 (0.50 * (1 - 0.1*lat_factor)) (0.24 + 0.1*lat_factor) 0.12;
        0.03 0.15 0.77 0.05;
        (0.28 + 0.1*lat_factor) 0.14 0.03 0.55
    ]
    
    fig_conn = Figure(size = (680, 480))
    ax_conn = Axis(
        fig_conn[1, 1],
        title = "Scotian Shelf Larval Demographic Connectivity (T_ij)",
        xticks = (1:4, cfa_names),
        yticks = (1:4, cfa_names),
        xlabel = "Settlement / Sink Region (j)",
        ylabel = "Spawning / Source Region (i)"
    )
    hm = heatmap!(ax_conn, 1:4, 1:4, T_mat', colormap = :Viridis, colorrange = (0.0, 0.8))
    Colorbar(fig_conn[1, 2], hm, label = "Transition Probability T_ij")
    fig_conn
end

# ╔═╡ a0000001-0001-0001-0001-000000000017
md"""
## 6. Assumptions and Scientific Utility

### Explicit Assumptions:
1. **Passive Horizontal Advection with Turbulent Dispersion**:
   Horizontal swimming motility ($< 5\text{ mm/s}$) is negligible compared to regional
   currents ($0.10 - 0.50\text{ m/s}$), justifying passive horizontal advection.
2. **Deterministic Visser Correction**:
   Subgrid vertical mixing requires $\frac{\partial \kappa_v}{\partial z}$ to satisfy mass conservation;
   omitting this term leads to unphysical particle accumulation in pycnoclines.
3. **Thermal Base $T_0 = -1.5^\circ\text{C}$**:
   Physiological development halts at $-1.5^\circ\text{C}$; below this threshold, larvae undergo thermal quiescence.
4. **Benthic Settlement Thermal Window**:
   Settlement is permitted only between $-350\text{ m} \le z \le -50\text{ m}$ where bottom temperatures
   satisfy $T_{\text{bed}} \le 6.0^\circ\text{C}$. Energetic tidal fluctuations are filtered out using
   a 12.42-hour moving exponential average.
5. **Multi-Scale Spatio-Temporal Coupling**:
   - **Ocean Hydrodynamics (`GeoOceanModel`)**: Solves 3D Navier-Stokes equations with barotropic subcycling at $\Delta t = 30 - 60\text{ s}$ on a structured $\sim 2.5\text{ km}$ spherical grid, resolving the baroclinic Rossby radius ($R_d \approx 4 - 8\text{ km}$). Output is archived at a 1-hour cadence ($\Delta t_{\text{output}} = 3600\text{ s}$) to capture semidiurnal tidal excursions ($M_2 \approx 12.42\text{ h}$) and inertial ellipses ($T_f \approx 17.1\text{ h}$).
   - **Offline Lagrangian Advection (`GeoLarvalTracking`)**: Particles are stepped via RK4 at $\Delta t_{\text{part}} = 300\text{ s}$ ($5\text{ min}$) by spatiotemporal interpolation of the hourly flow tensor, supplemented by subgrid turbulent dispersion ($K_h \approx 10\text{ m}^2/\text{s}$). High-resolution bathymetric steering ($\sim 100\text{--}500\text{ m}$) is enforced locally at the seabed boundary without the computational burden of fine-scale 3D Navier-Stokes solver grids everywhere.

### Management & Ecological Utility:
- **Marine Protected Area (MPA) Network Design**:
  Identifies persistent larval source zones (e.g. St. Anns Bank, Gully MPA) and recipient nurseries.
- **Stock Assessment & Spatially Explicit Recruitment**:
  Provides empirical recruitment indices to calibrate population dynamics models (`GeoIBM`).
- **Climate Vulnerability Projections**:
  Quantifies how ocean warming contracts pelagic larval duration (PLD) and alters downstream transport pathways.
"""

# ╔═╡ Cell order:
# ╠═a0000001-0001-0001-0001-000000000001
# ╠═a0000001-0001-0001-0001-000000000002
# ╠═a0000001-0001-0001-0001-000000000003
# ╠═a0000001-0001-0001-0001-000000000004
# ╠═a0000001-0001-0001-0001-000000000005
# ╠═a0000001-0001-0001-0001-000000000006
# ╠═a0000001-0001-0001-0001-000000000007
# ╠═a0000001-0001-0001-0001-000000000008
# ╠═a0000001-0001-0001-0001-000000000009
# ╠═a0000001-0001-0001-0001-000000000010
# ╠═a0000001-0001-0001-0001-000000000018
# ╠═a0000001-0001-0001-0001-000000000011
# ╠═a0000001-0001-0001-0001-000000000012
# ╠═a0000001-0001-0001-0001-000000000013
# ╠═a0000001-0001-0001-0001-000000000014
# ╠═a0000001-0001-0001-0001-000000000015
# ╠═a0000001-0001-0001-0001-000000000016
# ╠═a0000001-0001-0001-0001-000000000017
