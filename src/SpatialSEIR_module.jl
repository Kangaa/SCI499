##Spatial SEIR Source files
using DrWatson
@quickactivate "SCI499"

module SpatialSEIR

include(srcdir("MixingMatrices_src.jl"))
include(srcdir("JRB_Spatial_SEIR_src.jl"))

end