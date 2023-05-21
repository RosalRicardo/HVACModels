using  CairoMakie,DifferentialEquations, ModelingToolkit, Plots, GlobalSensitivity, Statistics

@variables t Tco(t)=15  

@parameters Cah(t)=4.5 fsw(t)=8.02e-5 ρw(t)=998 Cpw(t)=4.1868 Twi(t)=7 Two(t)=12 UAa(t)=0.04 To(t)=10 fsa(t)=0.192 ρa(t)=1.25 Cpa(t)=1.005 Tm(t)=22

D = Differential(t)

eqs = [
    D(Tco) ~ fsw*ρw*Cpw*(Twi-Two)+UAa*(To-Tco)+fsa*ρa*Cpa*(Tm-Tco)
]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,100.0)


ev_times = collect(0.0:1.0:100)
condition(u,t,integrator) = t ∈ ev_times
#affect!(integrator) = integrator.u[1] += 5*rand() print(integrator.p[7])
function affect!(integrator)
    integrator.p[1] += 6*rand()
    integrator.p[4] += 5*rand()
    println(integrator.p[4])
end

cb = DiscreteCallback(condition,affect!)

prob = ODEProblem(simpsys,[],tspan,callback=cb,tstops=ev_times)

sol = solve(prob)

Plots.plot(sol,title="Heating and Cooling coil model")
