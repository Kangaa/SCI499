function RK4(m::Model, x0::Array{Float64, 1}, t0::Float64, dt::Float64)::RK4 
    shape = x0|>length
    RK(
        k1 = zeros(shape),
        k2 = zeros(shape),
        k3 = zeros(shape),
        k4 = zeros(shape),
        k1_dt = zeros(shape),
        k2_dt = zeros(shape),
        k3_dt = zeros(shape),
        x1_k1 = zeros(shape),
        x1_k2 = zeros(shape),
        x1_k3 = zeros(shape),
        m = m,
        x = x0,
        t = t0,
        dt = dt
    )
end
