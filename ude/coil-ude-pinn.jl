using  CairoMakie,DifferentialEquations, ModelingToolkit, Plots, GlobalSensitivity, Statistics
using Flux, OrdinaryDiffEq, SciMLSensitivity, Plots

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

dados  = Array(sol)

# ===========================

# Multilayer FeedForward

ann = Chain(Flux.Dense(1,10,tanh), Dense(10,1))

p1,re = Flux.destructure(ann)
p2 = Float32[-2.0,1.1]
p3 = [p1;p_]
ps = Flux.params(p3)

function dudt_(du,u,p,t)
    Tco = u
    Cah, fsw, ρw, Cpw, Twi, Two, UAa, To, fsa, ρa, Cpa, Tm = p[end-12:end]
    NN = re(p[1:31])(u)[1]
    #du[1] = fsw*ρw*Cpw*(Twi-Two)+UAa*(To-Tco)+fsa*ρa*Cpa*(Tm-Tco)
    du[1] = fsw*ρw*Cpw*(Twi-Two)+UAa*(To-Tco[1])+NN
end

function dudt_ann(du,u,p,t)
    Tco = u
    Cah, fsw, ρw, Cpw, Twi, Two, UAa, To, fsa, ρa, Cpa, Tm = p[end-12:end]
    #du[1] = fsw*ρw*Cpw*(Twi-Two)+UAa*(To-Tco)+fsa*ρa*Cpa*(Tm-Tco)
    du[1] = re(p[1:31])(u)[1]
end


prob = ODEProblem(dudt_ann,u0,tspan,p3)
teste = concrete_solve(prob,Tsit5(),u0,p3,abstol=1e-8,reltol=1e-6,saveat=1)

teste.t

function predict_adjoint()
    Array(concrete_solve(prob,Tsit5(),u0,p3,saveat=1))
  end
  sum(abs2,dados-predict_adjoint())


  loss_adjoint() = sum(abs2,x-1 for x in predict_adjoint())
  loss_adjoint2() = sum(abs2,dados-predict_adjoint())
  loss_adjoint2()
  
  data = Iterators.repeated((), 500)
  opt = ADAM(0.05)
  iter = 0
  cb = function ()
    global iter += 1
    if iter % 50 == 0
      display(loss_adjoint2())
      display(Plots.plot(solve(remake(prob,p=p3,u0=u0),Tsit5(),saveat=0.1),ylim=(0,6)))
    end
  end

  Flux.train!(loss_adjoint2, ps, data, opt, cb = cb)

  Plots.plot(transpose(predict_adjoint()))
  Plots.plot!(transpose(dados))