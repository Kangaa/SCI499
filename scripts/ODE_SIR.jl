using DifferentialEquations
using Plots

function SIRModel!(du, u, p, t)
    S, I, R = u
    β, γ = p
    du[1] = -β*S*I
    du[2] = β*S*I-γ*I
    du[3] = γ*I
end 

p = [.1,.05]
u0 = [.99,.01,0]
tspan = (0.0, 200)

prob = ODEProblem(SIRModel!, u0, tspan, p)

sol = solve(prob)

plot(sol)