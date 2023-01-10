using  DifferentialEquations, ModelingToolkit, Plots

@variables t mᵣ(t)=0 W(t)=0

@parameters ηᵥ=0.98 ρ₁=998 Vₛ=0.98 Nrpm=4000 hs2=20 h1=20

D = Differential(t)

eqs = [
    D(mᵣ) ~ ηᵥ*ρ₁*Vₛ*(Nrpm/60)
    D(W) ~ mᵣ*((hs2-h1)/0.98)
]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,100.0)


ev_times = collect(0.0:1.0:100)

prob = ODEProblem(simpsys,[],tspan)

sol = solve(prob)

Plots.plot(sol)