function event_occurred(m::Model, event_ix::Int)
    ev_type::Int = event_ix  / m.numpatch
    ev_location::Int = event_ix  / m.numpatch
    if ev_type == 0 ## exposure
        m.s[ev_location] -= 1
        m.e[ev_location] += 1
    elseif ev_type == 1 ## Infection
        m.e[ev_location] -= 1
        m.i[ev_location] += 1
    elseif ev_type == 2 ## Recovery
        m.i[ev_location] -= 1
    end     
end