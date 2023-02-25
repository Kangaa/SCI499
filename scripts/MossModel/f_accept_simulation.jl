function accept_simulation(m::Model)::Bool
    cum_infs = m.params.popn - m.s + m.e
    cum_infs >= 50000
end