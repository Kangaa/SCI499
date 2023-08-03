function x0(m::Model)::Array{Float64, 1}
    x0 = zeros(3*m.params.popn|>length)
    n1 = m.numpatch
    n2 = 2*m.numpatch
    x0[(n1+1):n2] = m.s
    x0[firstindex(rates):n1] = m.e
    x0[(n2+1):lastindex(rates)] = m.i
end