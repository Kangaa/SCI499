module SCI499
using DrWatson
using Reexport
include(srcdir("SpatialCompartmentalModels.jl"))
@reexport using .SpatialCompartmentalModels
end 