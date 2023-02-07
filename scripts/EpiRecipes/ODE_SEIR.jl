using DifferentialEquations
using Plots

function seir_ode(dY, Y, p, t)
    S, E, I = Y
    β, σ, γ, μ = p

    dY[1]= μ - β*S*I - μ*S
    dY[2] = β*S*I - (σ + μ)*E
    dY[3] = σ*E - (γ+μ)*I
end



par=[520/365,1/60,1/30,774835/(65640000*365)]
init=[0.8,0.1,0.1]
tspan=(0.0,365.0)



prob = ODEProblem(seir_ode, init, tspan, par)

sol = solve(prob)

R=ones(1,size(sol,2))- sum(sol, dims = 1)

plot(sol.t, [sol', R'])