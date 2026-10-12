# ==============================================================================
#  store_reader.jl -- reading the realized simulation from its GeoData store
# ==============================================================================
#
#  Replaces `flow_reader.jl`, which read Oceananigans' JLD2 archive plus a separate grid
#  sidecar. That reader existed because the archive was the only thing a run produced; now the
#  run writes its slices into a Zarr store that carries its own coordinates, so a reader opens
#  one path and gets flow, axes, and times from it.
#
#  What the store provides (written by `GeoOceanModel.attach_output_callback!`):
#
#      hydrodynamics.zarr/
#        .zattrs          lon, lat, depth, dlon, dlat, fields, title
#        time.json        the time of each interval, once for the whole store
#        u/0  u/1  ...    one Zarr array per interval, shaped (n_lon, n_lat, n_depth)
#        T/0  T/1  ...
#
#  This package does not know how the store was produced. It knows the layout, and that is the
#  whole coupling between the two packages: `GeoOceanModel` writes it, this package reads it,
#  and neither imports the other. The layout is asserted in `test_store_reader.jl` rather than
#  assumed, so a change on either side shows as a failed read.
#
#  What is deliberately NOT here: checkpoints. Restart is the ocean half's concern and stays
#  in JLD2, the one format Oceananigans can checkpoint-restart from.
# ==============================================================================

using Zarr
using JSON3

"""
    _store_path(location) -> String

Resolve a location to the store directory.

A directory is the store itself; a path whose suffix is `.zarr` is used as given. Anything else
is an error naming what was received rather than a guess at what was meant, because a wrong
guess here reads a different simulation than the one asked for.
"""
function _store_path(location::AbstractString)
    p = abspath(location)
    isdir(p) && return p
    endswith(p, ".zarr") && isdir(p) && return p
    error("No hydrodynamic store at $(p). Expected a directory, or a `.zarr` store. This is " *
          "the location the run was handed, and a missing store cannot be substituted.")
end

"""
    _read_times(store) -> Vector{Float64}

The time of each written interval, in order.

Read from `time.json`, which the writer updates once per store rather than once per field. One
place holding the time axis is deliberate: two per-field copies could drift, and a field read
at its own recorded time but indexed by a shared time axis is a silent mismatch.
"""
function _read_times(store::AbstractString)
    tf = joinpath(store, "time.json")
    !isfile(tf) && error(
        "The store at $(store) has no `time.json`, so the times of its intervals are unknown. " *
        "A field with no time is not a simulation output.")
    isempty(tryparse(Float64, "0") ) && error("unreachable")
    parsed = JSON3.read(read(tf, String))
    out = Float64[]
    for v in parsed
        v isa Vector || (push!(out, Float64(v)); continue)
        for x in v
            push!(out, Float64(x))
        end
    end
    isempty(out) && error("`time.json` at $(store) records no times.")
    return out
end

"""
    _read_axes(store) -> (; lon, lat, depth, dlon, dlat)

The grid's interior cell centres, read from the store's own `.zattrs`.

Absent for a reason that is stated: these are what the sidecar used to carry, in a second file
next to the archive that nothing kept in step. They are now inside the store, and a store that
lacks them states which of the two is missing rather than defaulting to a spanned range --
defaulting would put every value in a track at the wrong place.
"""
function _read_axes(store::AbstractString)
    zattrs = joinpath(store, ".zattrs")
    !isfile(zattrs) && error(
        "The store at $(store) has no `.zattrs`, so its grid coordinates are unknown. A field " *
        "with no coordinates cannot be interpolated against.")
    a = JSON3.read(read(zattrs, String))
    get_vec(k) = [Float64(v) for v in a[k]]
    get_scalar(k) = Float64(a[k])
    has = all(haskey(a, k) for k in ("lon", "lat", "depth"))
    has || error(
        "The store at $(store) declares $(sort!(collect(keys(a)))) but not lon/lat/depth. " *
        "A store without its cell centres cannot be read.")
    return (lon = get_vec("lon"), lat = get_vec("lat"), depth = get_vec("depth"),
            dlon = get_scalar("dlon"), dlat = get_scalar("dlat"))
end

"""
    hydrodynamic_variables(store) -> Vector{String}

Which fields the store carries, from `.zattrs`. A reader that asked for a field the store does
not have would otherwise return zeros and call it a velocity.
"""
function hydrodynamic_variables(store::AbstractString)
    a = JSON3.read(read(joinpath(store, ".zattrs"), String))
    haskey(a, "fields") || return [n for n in readdir(store) if isdir(joinpath(store, n))]
    return [String(v) for v in a["fields"]]
end

"""
    open_hydrodynamics(location; fields = nothing) -> HydrodynamicStore

Open a realized simulation's store for reading.

`fields` defaults to everything the store carries, so a caller gets the whole run rather than
only what it happened to name. A request for a field the store does not have is an error
listing what it does have.

The result holds no flow arrays. The intervals are read on demand, which is what lets a
consumer hold a long run without holding it in memory.
"""
function open_hydrodynamics(location::AbstractString; fields = nothing)
    store = _store_path(location)
    times = _read_times(store)
    axes = _read_axes(store)
    available = hydrodynamic_variables(store)

    wanted = isnothing(fields) ? available :
             [String(f) for f in (fields isa Union{Tuple, AbstractVector} ? fields : (fields,))]
    unknown = [f for f in wanted if !(f in available)]
    isempty(unknown) || error(
        "The store at $(store) does not carry $(unknown). It carries $(available).")

    return HydrodynamicStore(store, available, times, axes)
end

"""
    HydrodynamicStore

A realized simulation's store, opened for reading: its location, its fields, its times, and
its axes. No flow data.

The separation is the point. A consumer is handed a location and reads the interval it needs,
rather than being handed every interval the run produced.
"""
struct HydrodynamicStore
    location::String
    variables::Vector{String}
    times::Vector{Float64}
    axes::@NamedTuple{lon::Vector{Float64}, lat::Vector{Float64},
                      depth::Vector{Float64}, dlon::Float64, dlat::Float64}
end

"""
    variable(store, name, m) -> Array{Float64, 3}

Interval `m` of one field, as `(n_lon, n_lat, n_depth)`.

Intervals are read from `store/<name>/<m>/`, each a Zarr array with its own `.zarray` metadata.
An out-of-range index is an error naming the interval count, not a clamped read: a track that
asked for interval 7 of a 4-interval run is a bug in the track, and returning the last interval
would hide it behind a plausible velocity.
"""
function variable(store::HydrodynamicStore, name::AbstractString, m::Integer)
    !(name in store.variables) && error(
        "The store at $(store.location) does not carry `$(name)`. It carries " *
        "$(store.variables).")
    (1 <= m <= length(store.times)) || error(
        "Interval $(m) is out of range: this run wrote $(length(store.times)) interval(s).")
    dir = joinpath(store.location, name, string(m - 1))
    !isdir(dir) && error(
        "Interval $(m) of `$(name)` is declared (time $(store.times[m])) but its directory " *
        "$(dir) is missing. The store and its time axis disagree.")
    return collect(Zarr.zopen(dir))
end

"""
    create_flow_interpolator(location; fields = ("u", "v", "w", "T"))

Build the `(lon, lat, z, t) -> (u, v, w, T)` function the larval model integrates against.

This replaces `create_flow_interpolator_from_jld2`. The difference is where the axes come from:
they are read from the store rather than re-derived from a sidecar, so a reader cannot disagree
with the run about where its cells are.

Bilinear in space and linear in time, with the nearest edge value outside the domain. A track
that wanders past the boundary still returns a real velocity rather than failing, which is
what `track_larval_cohort` expects of a flow field.
"""
function create_flow_interpolator(location::AbstractString;
                                  fields::Tuple = ("u", "v", "w", "T"))
    store = open_hydrodynamics(location; fields = fields)
    lon = store.axes.lon
    lat = store.axes.lat
    z = store.axes.depth
    t = store.times

    nx, ny, nz = length(lon), length(lat), length(z)
    nt = length(t)

    # Hold the arrays for the run. The store is read here rather than per step, because
    # `track_larval_cohort` samples it many times per interval, and re-reading a Zarr array
    # per sample turns the tracking step into I/O.
    data = Dict{String, Vector{Array{Float64, 3}}}()
    for name in fields
        data[name] = [variable(store, name, m) for m in 1:nt]
    end

    # The trilinear weights are the same for every field, so they are computed once per call
    # rather than once per field per step.
    _weight(vals, v) = begin
        i = something(findfirst(x -> x > v, vals), length(vals))
        i = clamp(i, 2, max(2, length(vals)))
        return (i, clamp((v - vals[i-1]) / (vals[i] - vals[i-1]), 0.0, 1.0))
    end
    _weight_single(vals, v) = begin
        (i, s) = _weight(vals, v)
        return (i, s)
    end
    _time_weight(v) = begin
        i = clamp(something(findfirst(x -> x > v, t), nt), 2, max(2, nt))
        θ = clamp((v - t[i-1]) / (t[i] - t[i-1]), 0.0, 1.0)
        return (i, θ)
    end

    return function flow_interpolator(lon_q::Real, lat_q::Real, z_q::Real, t_q::Real)
        i, sx = _weight(lon, Float64(lon_q))
        j, sy = _weight(lat, Float64(lat_q))
        k, sz = _weight(z, Float64(z_q))
        m, θ = _time_weight(Float64(t_q))

        function sample(name)
            arr = data[name]
            v0 = arr[m-1]; v1 = arr[m]
            a = 0.0
            for (ar, wf) in ((v0, 1.0 - θ), (v1, θ))
                a += wf * (
                    (1 - sx) * (1 - sy) * (1 - sz) * ar[i-1, j-1, k-1] +
                    sx       * (1 - sy) * (1 - sz) * ar[i,   j-1, k-1] +
                    (1 - sx) * sy       * (1 - sz) * ar[i-1, j,   k-1] +
                    sx       * sy       * (1 - sz) * ar[i,   j,   k-1] +
                    (1 - sx) * (1 - sy) * sz       * ar[i-1, j-1, k  ] +
                    sx       * (1 - sy) * sz       * ar[i,   j-1, k  ] +
                    (1 - sx) * sy       * sz       * ar[i-1, j,   k  ] +
                    sx       * sy       * sz       * ar[i,   j,   k  ])
            end
            return a
        end

        return (u = sample("u"), v = sample("v"), w = sample("w"), T = sample("T"))
    end
end

"""
    store_location(location::AbstractString) -> String

The location a consumer is handed. A path, and nothing else.
"""
store_location(location::AbstractString) = abspath(location)
