function step(rk::RK4)
    dt2 = rk.dt / 2 

    rk.k1 = f(rk.m, rk.x)
    rk.k1_dt = rk.k1*dt2
    rk.x_k1 = rk.k1*dt2

    rk.k2 = f(rk.m, rk.x_k1)
    rk.k2_dt = rk.k2*dt2
    rk.x_k2 = rk.x + rk.k2_dt

    rk.k3 = f(rk.m, rk.x_k2)
    rk.k3_dt = rk.k3*rk.dt
    rk.x_k3 = rk.x + rk.k3_dt
   
    rk.k4 = f(rk.m, rk.x_k3)
    rk.k2 *= 2.0
    rk.k3 *= 2.0
    rk.x += ((rk.k1+ rk.k2 + rk.k3 + rk.k4)* rk.dt / 6.0)
end