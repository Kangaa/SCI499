abstract type FracSelf end
struct Uniform{Float64} <: FracSelf end
struct Variable{Array} <: FracSelf end
