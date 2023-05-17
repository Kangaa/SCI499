# New source functions
module SpatialCompartmentalModels
using DrWatson
using Reexport

@reexport module MixingMatrices
    using DrWatson
    using DataFrames
    include(srcdir("MixingMatrices.jl"))
    export SpatialMixingMatrix
end

@reexport module CompartmentalModels
    using DrWatson
    using DataFrames
    include(srcdir("CompartmentalModels.jl"))
    export simulate!
end

end