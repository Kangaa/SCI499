
module SCI499
using DrWatson

include(srcdir("MixingMatrices_src.jl"))
include(srcdir("JRB_Spatial_SEIR_src.jl"))

using Reexport
@reexport using .SpatialSEIR, .MixingMatrices, Shapefile, Tables, DataFrames
end