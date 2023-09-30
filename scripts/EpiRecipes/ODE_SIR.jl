using DifferentialEquations
using Plots

function SIRModel!(du, u, p, t)
    S, I, R = u
    β, γ = p
    du[1] = -β*S*I * fill(0.5,4,4)
    du[2] = β*S*I-γ*I
    du[3] = γ*I
end 

p = [.1,.05]
u0 = [fill(.99,4), fill(.01, 4)]
tspan = (0.0, 200)

prob = ODEProblem(SIRModel!, u0, tspan, p)

sol = solve(prob)

plot(sol)