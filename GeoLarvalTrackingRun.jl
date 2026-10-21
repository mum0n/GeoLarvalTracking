# ==============================================================================
#  GeoLarvalTrackingRun.jl -- CLI for the Lagrangian half
# ==============================================================================
#
#  This file holds CLI interaction only: environment bootstrap, flag parsing, and dispatch.
#  Everything that defines behaviour -- the bathymetry/domain derivation, the segments, the
#  figure families -- is in `src/`, because a script that defines the behaviour is a script
#  nobody can test, and `Pkg.test` cannot reach a function a script defines.
#
#  The verbs, deliberately few. A parameter a run should repeat belongs in the TOML, not the
#  command line; a flag for it is a flag whose value nobody records:
#
#    rerun   Run the tracking again, reusing the existing store rather than rebuilding it.
#    force   Re-run and overwrite the store, so a repeat of an existing run replaces it.
#    quick   A short unsaved run for checking that the configuration still works.
#    viz     Render the figures from the stored results.
#
#  On the fourth ("is there a better way?"): `Pkg.test` is the better route for *checking*, and
#  it needs a `test/` suite first, which is the same prerequisite the driver move needed. The
#  CLI here is for producing and re-reading results, not for asserting anything about them.
#
#  Usage:
#      julia --project=. GeoLarvalTrackingRun.jl --help
#      julia --project=. GeoLarvalTrackingRun.jl --rerun
# ==============================================================================

# Captured at top level: `@__FILE__` is only valid at file scope.
const _THIS_SCRIPT = @__FILE__

# Auto-discover and activate this package's environment when run as a script. Guarded so the
# bootstrap is skipped when this file is `include`d from a test, and wrapped in a `let` so the
# search's assignments have a local scope instead of ambiguous soft-scoped script globals.
if !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == _THIS_SCRIPT &&
   Base.find_package("GeoLarvalTracking") === nothing
    import Pkg
    let curr = @__DIR__
        while true
            proj = joinpath(curr, "Project.toml")
            if isfile(proj) && occursin("name = \"GeoLarvalTracking\"", read(proj, String))
                Pkg.activate(curr)
                break
            end
            parent = dirname(curr)
            parent == curr && break
            curr = parent
        end
    end
end

using GeoLarvalTracking
using GeoLarvalTracking: configuration_to_options, options_to_configuration,
    save_configuration, load_configuration, resolve_config_name, run_segment_tracking,
    run_segment_metrics, run_segment_visualize, render_all_figures, HydrodynamicOptions
using GeoData

"""
    display_help(io::IO = stdout)

The four verbs and the two flags that modify them.
"""
function display_help(io::IO = stdout)
    println(io, """
    GeoLarvalTrackingRun.jl -- larval tracking through a realized ocean simulation.

    Usage:
      julia --project=. GeoLarvalTrackingRun.jl [flags]

    Verbs:
      --rerun      Run the tracking again, reusing the existing hydrodynamic store.
      --force      Same, and overwrite the stored results rather than skipping them.
      --overwrite  Alias for --force: overwrite stored results and re-track from t = 0.
      --quick      A short unsaved run to check the configuration still works. Writes
                   nothing that a later --rerun would need to skip.
      --viz        Render the figures from the stored results, without re-tracking.

    Clearing what a previous run left behind:
      --blank_state  Wipe ALL of this run's state, the downloaded store included: the
                    trajectory store, the trajectory archives, the restart checkpoints, and
                    the ingested inputs. The next run fetches the inputs again. Asks
                    before it deletes.
      --force and --overwrite also ask before they delete, and they delete the MODEL
      RESULTS only: the trajectory store, the trajectory archives and the restart
      checkpoints. The downloaded store is kept, so the run does not fetch it again.

    Configuration:
      --config=FILE        A path, or a name resolved against configs/.
      --hydro-file=PATH    The realized ocean simulation to track through. REQUIRED for
                          --rerun, --force and --quick: this is the only coupling with
                          GeoOceanModel, and it is a path. Its metadata states the domain, so
                          the ocean half's configuration is never read here.
      --save-config[=D]    Write the fully-resolved configuration, then exit.

    Two configurations:
      default       Cape Breton  demonstration -- small cohort, generic placeholder
                    biology; run against the keyless ocean default.toml store.
      snowcrab      Snow Crab production -- 500-particle cohort over the 60-day pelagic
                    larval duration with snow-crab biology, run against the ocean
                    scotian_shelf.toml store (GeoOceanModel/configs/scotian_shelf.toml).

    There are no per-run override flags. A particle count or a duration belongs in the TOML so
    the value is recorded with the run; a flag for it is a flag whose value nobody writes down,
    and a re-run from the config would not reproduce what the flag did.
    """)
    return nothing
end

"""
    parse_args(args::Vector{String}) -> Dict{String, Any}

The whole flag vocabulary in one place. An unrecognised flag is an error rather than silence, so
a typo cannot quietly select the default verb.

Exactly one verb is accepted. A command line naming two is a caller who has not decided, and
guessing which they meant is how a long run happens instead of a figure.
"""
function parse_args(args::Vector{String})
    out = Dict{String, Any}(
        "verb" => nothing, "config" => "", "hydro_file" => "", "save_config" => nothing,
        "blank_state" => false, "help" => false)
    verbs = Set(["--rerun", "--force", "--overwrite", "--quick", "--viz"])
    for a in args
        if a in verbs
            verb_name = a == "--overwrite" ? "force" : replace(a, "--" => "")
            out["verb"] !== nothing && error(
                "Two verbs given: --$(verb_name) and --$(out["verb"]). Name one; " *
                "which one was meant is not a thing this script should guess.")
            out["verb"] = verb_name
        elseif a in ("--blank_state", "--blank-state")
            out["blank_state"] = true
        elseif a in ("--help", "-h")
            out["help"] = true
        elseif a == "--save-config"
            out["save_config"] = ""
        elseif startswith(a, "--save-config=")
            out["save_config"] = split(a, "=", limit = 2)[2]
        elseif startswith(a, "--config=")
            out["config"] = split(a, "=", limit = 2)[2]
        elseif startswith(a, "--hydro-file=")
            out["hydro_file"] = split(a, "=", limit = 2)[2]
        else
            error("Unknown flag: $(a). Run with --help.")
        end
    end
    return out
end

"""
    resolve_config_file(flag::AbstractString) -> String

Which configuration this run uses, as a *path*. An empty flag is this half's own
`configs/default.toml`; anything else is resolved through `resolve_config_path`.
"""
function resolve_config_file(flag::AbstractString)
    isempty(strip(flag)) && return joinpath(@__DIR__, "configs", "default.toml")
    return resolve_config_path(String(flag))
end

"""
    resolve_config_path(flag::AbstractString) -> String

A file path is used directly; a bare name is looked up under this package's `configs/`
directory. A name that resolves to nothing is an error naming what was tried, rather than a
fall-back to a default the operator did not choose.
"""
function resolve_config_path(flag::AbstractString)::String
    f = String(flag)
    isfile(f) && return abspath(f)
    name = endswith(f, ".toml") ? f : f * ".toml"
    named = joinpath(@__DIR__, "configs", name)
    isfile(named) && return named
    error("Configuration file not found: $(f). Pass a path, or a name under configs/ " *
          "(for example: default, snowcrab).")
end

"""
    clear_run_state(opts::HydrodynamicOptions; blank_state::Bool, input = stdin) -> Bool

Clear what a previous run of this configuration left behind, after warning about it and
being told to go ahead. Returns `false` when the operator declined, in which case the
caller must not run.

Two kinds of state live in this run's output directory, and the two flags mean different
things because of it:

  * the MODEL RESULTS -- the trajectory store, the trajectory archives and the restart
    checkpoints. `--force` clears these and keeps the rest.
  * the DOWNLOADED STORE -- this run's `inputs/`: the coastal polygons and the vertical
    grid. `--blank_state` clears these as well, so the next run fetches them again.

The shared GeoData store under the lakehouse is never touched: it holds the ocean half's
flow solution, which this half tracks through and does not own.

`input` is the stream the answer is read from: `stdin` in a run, and a buffer when the
run is being driven from a script.
"""
function clear_run_state(opts::HydrodynamicOptions; blank_state::Bool = false,
                        input::IO = stdin)
    out_dir = opts.output_dir
    tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
    results = [opts.zarr_path,
               joinpath(out_dir, "larval_trajectories.jld2"),
               joinpath(out_dir, "larval_trajectories$(tag).jld2")]
    result_globs = isempty(strip(opts.checkpoint_dir)) ? String[] :
                   [joinpath(opts.checkpoint_dir, "*.jld2")]
    inputs = [opts.input_dir]
    return GeoData.purge_run_state(out_dir; results, result_globs, inputs,
                                   blank_state = blank_state, input = input)
end

function main(args::Vector{String})
    flags = parse_args(args)
    flags["help"] && return display_help()

    cfg_file = resolve_config_file(flags["config"])
    # "Overwrite" is a property of the RUN, and it is carried in the options rather than set
    # on them: `HydrodynamicOptions` is immutable, and a resolved options struct that could
    # be edited after the fact would not be a resolved struct. The configuration's own
    # `overwrite` is left alone when the flags do not ask for one.
    overwrite_requested = flags["verb"] == "force" || flags["blank_state"]
    opts = try
        # The configuration's own name is passed on: it is what names this run's directory
        # in the lakehouse when the configuration does not name one itself.
        configuration_to_options(load_configuration(cfg_file);
                                 config_name = resolve_config_name(cfg_file),
                                 (overwrite_requested ? (; overwrite = true) : (;))...)
    catch err
        println("Configuration error: $(sprint(showerror, err))")
        return display_help()
    end

    # --save-config works without a verb, as it does in the ocean half: resolve and save the
    # fully-resolved configuration, then exit. (A run needs a verb; inspecting the config does
    # not.)
    if !isnothing(flags["save_config"])
        if flags["blank_state"]
            error("--save-config writes the resolved configuration and exits; it cannot also " *
                  "clear a previous run's state. Run them one at a time.")
        end
        dest = isempty(flags["save_config"]) ? cfg_file : flags["save_config"]
        save_configuration(options_to_configuration(opts), dest)
        println("Resolved configuration written to $(dest)")
        return nothing
    end

    isnothing(flags["verb"]) && return display_help()

    if flags["blank_state"] && flags["verb"] == "viz"
        error("--viz renders the figures from the stored results, so it cannot also wipe " *
              "them. Use --rerun or --force with --blank_state.")
    end
    if flags["blank_state"] && flags["verb"] == "quick"
        error("--quick writes nothing that a later --rerun would need, so it cannot also " *
              "wipe a previous run's state. Use --rerun or --force with --blank_state.")
    end

    println("="^73)
    println(" GeoLarvalTracking")
    println("="^73)
    println(" Verb:       --$(flags["verb"])")
    println(" Config:     $(cfg_file)")
    println(" Hydro file: $(isempty(flags["hydro_file"]) ? "(not given)" : flags["hydro_file"])")
    println("="^73)

    # The one coupling with the ocean half, and it is a path. Everything about the ocean --
    # domain, grid, bathymetry -- comes from that file's metadata, so the ocean half's own
    # configuration is never opened here.
    hydro = flags["hydro_file"]
    needs_hydro = flags["verb"] in ("rerun", "force", "quick")
    if needs_hydro && isempty(strip(hydro))
        error("--$(flags["verb"]) needs the location of a realized ocean simulation. Pass " *
              "--hydro-file=PATH, produced by GeoOceanModelRun.jl.")
    end

    if flags["verb"] == "rerun" || flags["verb"] == "force" || flags["verb"] == "quick"
        # What a previous run left behind is cleared BEFORE the tracking starts, so the run
        # that follows sees the state it is meant to see. --force clears the model results
        # and keeps the downloaded store; --blank_state clears both, so the inputs are
        # fetched again. Either way the operator is asked first.
        if flags["verb"] == "force" || flags["blank_state"]
            cleared = clear_run_state(opts; blank_state = flags["blank_state"])
            if !cleared
                println("Aborted: the existing state was left untouched.")
                return nothing
            end
        end
        # A rerun reads what the existing store says and skips what it already covers; --force /
        # --overwrite replaces stored results and re-tracks from t = 0. --quick runs tracking
        # without metrics, for checking that the configuration still works.
        trajs = run_segment_tracking(opts = opts, hydro_file = hydro)
        if flags["verb"] != "quick"
            run_segment_metrics(opts = opts, trajectories = trajs)
        end
    elseif flags["verb"] == "viz"
        render_all_figures(opts = opts)
    end
    return nothing
end

is_cli = !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == _THIS_SCRIPT
is_cli && main(ARGS)
