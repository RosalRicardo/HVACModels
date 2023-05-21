using Flux, Plots, Statistics

using  CairoMakie,DifferentialEquations, ModelingToolkit, Plots, GlobalSensitivity, Statistics

@variables t Tco(t)=15  

u0 = Float64[15.0]

@parameters Cah(t)=4.5 fsw(t)=8.02e-5 ρw(t)=998 Cpw(t)=4.1868 Twi(t)=7 Two(t)=12 UAa(t)=0.04 To(t)=10 fsa(t)=0.192 ρa(t)=1.25 Cpa(t)=1.005 Tm(t)=22

p_ = [4.5,8.02e-5,998,4.1868,7,12,0.04,10,0.192,1.25,1.005,22]

D = Differential(t)

eqs = [
    D(Tco) ~ fsw*ρw*Cpw*(Twi-Two)+UAa*(To-Tco)+fsa*ρa*Cpa*(Tm-Tco)
]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,200.0)


ev_times = collect(0.0:1.0:100)
condition(u,t,integrator) = t ∈ ev_times
#affect!(integrator) = integrator.u[1] += 5*rand() print(integrator.p[7])
function affect!(integrator)
    integrator.p[1] += 3*rand()
    integrator.p[4] += 3*rand()
    println(integrator.p[4])
end

cb = DiscreteCallback(condition,affect!)

prob = ODEProblem(simpsys,[],tspan)

sol = solve(prob,saveat=1)

Plots.plot(sol,title="Heating and Cooling coil model")

dados = Array(sol)



# ====

NNODE = Chain(x -> [x],
           Dense(1,32,tanh),
           Dense(32,1),
           first)

NNODE(1.0)

g(t) = t*NNODE(t) + 15

ϵ = sqrt(eps(Float32))

loss() = mean(abs2(((g(t+ϵ)-g(t))/ϵ) - cos(2π*t)) for t in 0:1f-2:1f0)

loss2() = mean(abs2(((g(t+1)-g(t))/1) - dados[t]) for t in 1:1:200)

loss_det() = mean(abs2(((g(t+1)-g(t))/ϵ) - 8.02e-5*998*4.1868*(7-12)+0.04*(10-dados[t])+0.192*1.25*1.005*(22-dados[t])) for t in 1:1:200)

mean(abs2(((g(t+ϵ)-g(t))/ϵ) - dados[t]) for t in 1:1:200)

g()

opt = Flux.Descent(0.01)
data = Iterators.repeated((), 5000)
iter = 0
cb = function () #callback function to observe training
  global iter += 1
  if iter % 500 == 0
    display(loss2())
  end
end
display(loss2())

Flux.train!(loss2, Flux.params(NNODE), data, opt; cb=cb)

t = 1:1:200
Plots.plot(t,g.(t),label="NN")

g(1+)
Plots.plot!(t,transpose(dados)[1:200], label = "True")

