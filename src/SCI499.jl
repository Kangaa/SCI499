
module SCI499


using DrWatson

include(srcdir("SpatialSEIR_src.jl"))

using Reexport
@reexport using .SpatialSEIR, .MixingMatrices, Shapefile, Tables, DataFrames, CSV
end