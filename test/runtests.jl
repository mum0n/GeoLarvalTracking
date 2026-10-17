# ==============================================================================
# GeoLarvalTracking Comprehensive Test Suite Orchestrator
# ==============================================================================

import Pkg

const RUNNING_UNDER_PKG_TEST = basename(dirname(Base.active_project())) == "test"
if !RUNNING_UNDER_PKG_TEST
    Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
    Pkg.instantiate(; io = devnull)
end

# Usage:
#   julia --project=. test/runtests.jl                # Run every segment
#   julia --project=. test/runtests.jl store          # Run only store segment
#   julia --project=. test/runtests.jl driver         # Run only driver segment
#
# Under `Pkg.test()` the whole suite runs; pass segments as test_args:
#   Pkg.test(test_args = ["store"])
# ==============================================================================

using Test
using GeoLarvalTracking
using GeoData
using JSON3

const AVAILABLE_SEGMENTS = Dict{Symbol, @NamedTuple{file::String, desc::String}}(
    :store => (
        file = "test_store.jl",
        desc = "Hydrodynamics Store Reader, Interpolator & Verification"
    ),
    :driver => (
        file = "test_driver.jl",
        desc = "Pipeline Driver Segments & Surface API Export"
    )
)

"""
    parse_requested_segments(args::Vector{String}) -> Vector{Symbol}

Resolves CLI segment arguments to segment keys.
"""
function parse_requested_segments(args::Vector{String})::Vector{Symbol}
    if isempty(args)
        return collect(keys(AVAILABLE_SEGMENTS))
    end

    selected = Symbol[]
    for a in args
        key = Symbol(lowercase(strip(a)))
        if haskey(AVAILABLE_SEGMENTS, key)
            push!(selected, key)
        else
            avail = sort(string.(collect(keys(AVAILABLE_SEGMENTS))))
            @warn "Unrecognized test segment: '$a'. Available segments: $(join(avail, ", "))"
        end
    end

    return isempty(selected) ? sort!(collect(keys(AVAILABLE_SEGMENTS))) : unique(selected)
end

requested_segments = parse_requested_segments(ARGS)

@testset "GeoLarvalTracking Comprehensive Suite" begin
    for seg_key in requested_segments
        info = AVAILABLE_SEGMENTS[seg_key]
        println("\n" * "="^80)
        println(">>> Running Segment: $(info.desc) ($(info.file))")
        println("="^80)
        include(joinpath(@__DIR__, info.file))
    end
end
