"""
    visualization.jl

The output families for the Lagrangian half: the larval dispersal figures and the shared
figure helpers both halves use.

This is the larval package's counterpart to `GeoOceanModel`'s `output/visualization.jl`, which
carries the ocean figure families instead. `shared_figures.jl` is here because it is the
smallest family that both packages use, and this package is the one that needs it under the
name it is used under.
"""

using CairoMakie
import CairoMakie: record

include(joinpath(@__DIR__, "figures", "larval_disperal.jl"))
include(joinpath(@__DIR__, "figures", "shared_figures.jl"))
