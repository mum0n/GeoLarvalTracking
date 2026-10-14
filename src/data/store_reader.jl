# ==============================================================================
#  store_reader.jl -- reading the realized simulation from its GeoData store
# ==============================================================================
#
#  Replaces `flow_reader.jl`, which read Oceananigans' JLD2 archive plus a separate grid
#  sidecar. That reader existed because the archive was the only thing a run produced; now the
#  run writes its slices into a GeoData store that carries its own coordinates, so a reader
#  opens one path and gets flow, axes, and times from it.
#
#  What the store provides (written by `GeoOceanModel.attach_output_callback!`):
#
#      hydrodynamics.zarr/          a GeoDataset
#        .zattrs                    title, created_utc, and the seabed the run immersed
#        lon, lat, depth, time      coordinate arrays, in the dataset's own `coords`
#        T, u, v, w                 fields, each (lon, lat, depth, time)
#
#  Reading it is `GeoData.geoload(store)`. That is the whole point of the store being a GeoData
#  store: an earlier writer put each interval in its own directory (`T/0`, `T/1`), which GeoData's
#  Zarr reader cannot open -- it discovers variables as the root group's child arrays -- so every
#  consumer walked Zarr itself. This did. Now it does not.
#
#  This package does not know how the store was produced. It knows the layout, and that is the
#  whole coupling between the two packages: `GeoOceanModel` writes it, this package reads it,
#  and neither imports the other. The layout is asserted in `test_store_reader.jl` rather than
#  assumed, so a change on either side shows as a failed read.
#
#  What is deliberately NOT here: checkpoints. Restart is the ocean half's concern and stays
#  in JLD2, the one format Oceananigans can checkpoint-restart from.
# ==============================================================================

using GeoData
using GeoData: GeoDataset

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

The time of each written interval, in order, from the store's own `time` coordinate.

One axis for the whole store is deliberate: a per-field copy could drift, and a field read at
its own recorded time but indexed by a shared time axis is a silent mismatch. The coordinate
is read here from a `geoload`ed dataset, so it is the same axis the ocean half wrote.
"""
function _read_times(store::AbstractString)
    ds = _geoload(store)
    haskey(ds.coords, :time) || error(
        "The store at $(store) has no `time` coordinate, so the times of its intervals are " *
        "unknown. A field with no time is not a simulation output.")
    out = collect(Float64, ds.coords[:time].data)
    isempty(out) && error("The store at $(store) records no times.")
    return out
end

"""
    _read_axes(store) -> (; lon, lat, depth)

The grid's interior cell centres, from the store's own dimension records.

Present for the reason the sidecar used to be: these are what a reader needs to place a value
on Earth, and a store that does not state cannot be interpolated against. They come from the
dataset rather than from a second file that nothing kept in step, and a store that declares
fewer than three axes errors naming which, rather than defaulting to a spanned range --
defaulting would put every value in a track at the wrong place.
"""
function _read_axes(store::AbstractString)
    ds = _geoload(store)
    out = (lon = nothing, lat = nothing, depth = nothing)
    for (k, name) in ((:lon, "longitude"), (:lat, "latitude"), (:depth, "depth"))
        haskey(ds.dims, k) || haskey(ds.dims, Symbol(name)) || error(
            "The store at $(store) declares $(sort!(collect(keys(ds.dims)))) but not $(k). " *
            "A store without its cell centres cannot be read.")
        c = haskey(ds.dims, k) ? collect(Float64, ds.dims[k].coords) :
                                 collect(Float64, ds.dims[Symbol(name)].coords)
        out = merge(out, NamedTuple{(k,)}((c,)))
    end
    return out
end

"""
    hydrodynamic_variables(store) -> Vector{String}

Which fields the store carries. A reader that asked for a field the store does not have would
otherwise return zeros and call it a velocity.
"""
function hydrodynamic_variables(store::AbstractString)
    return sort!(collect(keys(_geoload(store).variables)))
end

"""
    _geoload(store) -> GeoDataset

The store as a GeoData dataset, opened through GeoData.

An empty variable set is an error here rather than at the first read: it is what a store that
is not GeoData's shape produces -- `geoload` succeeds and finds nothing -- and the failure used
to surface three call sites later as "does not carry T", naming the wrong layer.
"""
function _geoload(store::AbstractString)
    ds = GeoData.geoload(_store_path(store); backend = :zarr)
    isempty(ds.variables) && error(
        "The store at $(store) opens as a GeoData dataset with no variables. Either it is " *
        "empty, or it was written by an older writer that put each interval in its own " *
        "directory (`T/0`, `T/1`), which GeoData cannot read. Re-run the simulation to " *
        "produce a GeoData store.")
    return ds
end

"""
    open_hydrodynamics(location; fields = nothing) -> HydrodynamicStore

Open a realized simulation's store for reading.

`fields` defaults to everything the store carries, so a caller gets the whole run rather than
only what it happened to name. A request for a field the store does not have is an error
listing what it does have.

The result holds the opened dataset. The intervals are read on demand, which is what lets a
consumer hold a long run without holding it in memory.
"""
function open_hydrodynamics(location::AbstractString; fields = nothing)
    store = _store_path(location)
    ds = _geoload(store)
    times = collect(Float64, ds.coords[:time].data)
    axes = _read_axes(store)
    available = sort!(collect(keys(ds.variables)))

    wanted = isnothing(fields) ? available :
             [String(f) for f in (fields isa Union{Tuple, AbstractVector} ? fields : (fields,))]
    unknown = [f for f in wanted if !(f in available)]
    isempty(unknown) || error(
        "The store at $(store) does not carry $(unknown). It carries $(available).")

    return HydrodynamicStore(store, available, times, axes, ds)
end

"""
    HydrodynamicStore

A realized simulation's store, opened for reading: its location, its fields, its times, its
axes, and the dataset they came from. No flow data.

The separation is the point. A consumer is handed a location and reads the interval it needs,
rather than being handed every interval the run produced.
"""
struct HydrodynamicStore
    location::String
    variables::Vector{String}
    times::Vector{Float64}
    axes::@NamedTuple{lon::Vector{Float64}, lat::Vector{Float64},
                      depth::Vector{Float64}}
    dataset::GeoDataset
end

"""
    variable(store, name, m) -> Array{Float64, 3}

Interval `m` of one field, as `(n_lon, n_lat, n_depth)`.

Read from the dataset's own array, at `[:, :, :, m]`. An out-of-range index is an error naming
the interval count, not a clamped read: a track that asked for interval 7 of a 4-interval run
is a bug in the track, and returning the last interval would hide it behind a plausible
velocity.
"""
function variable(store::HydrodynamicStore, name::AbstractString, m::Integer)
    !(name in store.variables) && error(
        "The store at $(store.location) does not carry `$(name)`. It carries " *
        "$(store.variables).")
    (1 <= m <= length(store.times)) || error(
        "Interval $(m) is out of range: this run wrote $(length(store.times)) interval(s).")
    v = store.dataset.variables[name]
    size(v.data, 4) >= m || error(
        "Interval $(m) is declared by the time axis but field `$(name)` has only " *
        "$(size(v.data, 4)) interval(s). The store and its time axis disagree.")
    # Copied to the host rather than returned as a disk-backed view: the interpolator reads
    # this slice many times per tracking step, and a view is a decompression each time.
    return Array{Float64}(v.data[:, :, :, m])
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
