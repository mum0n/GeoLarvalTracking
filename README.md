# GeoLarvalTracking

Snow-crab larval-connectivity modelling. It tracks snow-crab larvae through a realized ocean
simulation and derives recruitment/connectivity metrics.

## What it does

`GeoLarvalTrackingRun.jl` dispatches (`src/driver.jl`, library code) to:

- **`run_segment_tracking`** — release a larval cohort into an ocean store, advect/disperse/DVM/
  degree-day-molt the particles through it, and store trajectories.
- **`run_segment_metrics`** — recruitment/connectivity from the trajectories.
- **`render_all_figures`** — figures.

### What it reads

The advective field (`src/data/store_reader.jl`, `create_flow_interpolator`) opens a
hydrodynamic store and reads `u,v,w,T` per interval with the store's own cell-centre
coordinates. The spawn domain is taken from the store's own axes when present, CFA boundaries
otherwise.

That store is **input data produced elsewhere**: a `hydrodynamics.zarr` GeoData store written by
an ocean run — [`GeoOceanModel`](../GeoOceanModel/README.md) produces the ones this package's
configs default to. It is read-only here. The store's own coordinates are used rather than a
grid declared in this package, so the two cannot drift apart.

## Running

```
julia --project=. GeoLarvalTrackingRun.jl --rerun \
  --hydro-file=<ocean store> \
  --config=default            # or --config=snowcrab
```

`--rerun` **or** `--force` run tracking; without a `--hydro-file` it needs one — this is the one
coupling, and it is a path, not a shared config. `--verbose`/`--viz` render figures from the
stored results.

## Configurations

Two, in `configs/` (the biology schema sections — `biology`, `dvm`,
`molting_and_settlement` — are the larval model's own):

| config | purpose |
|---|---|
| `default.toml` | Cape Breton demonstration — small cohort, generic placeholder biology, tracked through the keyless ocean `default.toml` store. |
| `snowcrab.toml` | Snow Crab production — 500-particle cohort over the 60-day pelagic larval duration, with snow-crab biology/DVM/degree-day parameters. |

- **`[storage]`** is where the run reads and writes: `output_dir` (this package's results
  directory) and `ocean_model` (the realized ocean store to track through — the
  `--hydro-file` default). Both default into `lakehouse_dir("Geo{Larval,Ocean}Model", …)`;
  a command-line `--hydro-file` overrides. `ocean_model` names the GeoData store itself —
  the `hydrodynamics.zarr` directory the ocean half wrote — not the run directory holding
  it.
- A configuration that names no output directory of its own still writes into the
  lakehouse, under its own name. A run never lands in the package root: that is outside
  the lakehouse, where the other half never looks and no later run finds it.
- No `--config` ⇒ `configs/default.toml`.

### Shared Lakehouse Configuration

`GeoLarvalTracking` resolves all lakehouse directories via `GeoData.lakehouse_dir(...)`.
Precedence hierarchy:
* `ENV["GEODATA_LAKEHOUSE"]` (global override).
* `LocalPreferences.toml` (`[GeoData] lakehouse_dir`).
* `~/.geodata/storage.toml` (`[storage] store_root`).
* Proximity fallback (`projects/lakehouse/`).

To configure the lakehouse persistently:
```julia
using GeoData
GeoData.set_lakehouse_dir!("/path/to/shared/lakehouse")
```

## Results

Trajectories, metrics and figures are written to this package's results store
(`[storage] output_dir`, a Zarr store via the analysis layer) — separate from the input
ocean store, which is read only.

## Notes / known gaps

- `run_segment_tracking(opts=…)` takes the ocean model's options type; the biology configs carry
  the larval biology keys, and the ocean configuration is never loaded here — the store's own
  metadata states the domain.
- `Build` (`src/driver.jl`) is a deliberate marker: a local helper of the umbrella driver
  whose body was never ported. It is uncalled; a run that needs it should port it explicitly.
- `SimulationResult` (the ocean's `open_simulation_result`) is a separate metadata handover
  that is defined but not yet written by a run; the live coupling is the store path.

## Documentation

- [`docs/GeoLarvalTracking.md`](docs/GeoLarvalTracking.md) — the full reference manual
  (configuration, biology, analysis, storage, visualization).
- [`docs/snow_crab_larval_connectivity_paper.md`](docs/snow_crab_larval_connectivity_paper.md)
  — the snow-crab larval connectivity paper.
