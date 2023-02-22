using DiffEqFlux, Optimization, OptimizationFlux, DifferentialEquations, LinearAlgebra
k, α, β, γ = 1, 0.1, 0.2, 0.3
tspan = (0.0,10.0)

function dxdt_train(du,u,p,t)
  du[1] = u[2]
  du[2] = -k*u[1] - α*u[1]^3 - β*u[2] - γ*u[2]^3
end

u0 = [1.0,0.0]
ts = collect(0.0:0.1:tspan[2])
prob_train = ODEProblem{true}(dxdt_train,u0,tspan)
data_train = Array(solve(prob_train,Tsit5(),saveat=ts))

A = [LegendreBasis(10), LegendreBasis(10)]
nn = TensorLayer(A, 1)

f = x -> min(30one(x),x)

function dxdt_pred(du,u,p,t)
  du[1] = u[2]
  du[2] = -p[1]*u[1] - p[2]*u[2] + f(nn(u,p[3:end])[1])
end

α = zeros(102)

prob_pred = ODEProblem{true}(dxdt_pred,u0,tspan)

