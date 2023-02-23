using  CairoMakie,DifferentialEquations, ModelingToolkit, Plots, GlobalSensitivity, Statistics

# variables

# Cz  - Overall thermal capacitance of the zone - 47.1 kJ/C
# Tz  - Zone Temperature - C
# Fsa - volume flow rate of supply fan - 0.192 m3/s
# ρa  - density of air - 1.25 kg/m3 
# Cpa - specific heat of supply air - 1.005 kJ/kgC
# Tsa - temperature of supply air - C
# Uw1 - overall heat transfer coefficient (East and west walls) - 2 w/m2C
# Uw2 - overall heat transfer coefficient (North and South walls) - 2 w/m2C
# Ur  - overall heat transfer coefficient (roof) - 1 w/m2C
# Aw1 - Area (East and West walls) - m2
# Aw2 - Area (North and South walls) - m2 
# Ar  - Area (Roof) - m2
# Tw1 - temperature (East and West walls) - C
# Tw2 - temperature (North and South walls) - C
# Tr  - temperature (Roof) - C
# q   - heat gain for occupants and lights - W

@variables t Tz(t)=40 Tw1(t)=20 Tw2(t)=20 Tr(t)=25 Wz(t)=0.5

@parameters Cz=47.1e3 Fsa=0.192  ρa=1.25 Cpa=1.005 Tsa=8 Uw1=2 Uw2=2 Ur=1 Aw1=9 Aw2=12 Ar=9 q=3000 To=21 Cw1=70 Cw2=60 Cr=80 Vz=36 Ws=0.02744 P=0.08

D = Differential(t)

eqs = [D(Tz) ~ (Fsa*ρa*Cpa*(Tsa-Tz)+2*Uw1*Aw1*(Tw1-Tz)+Ur*Ar*(Tr-Tz)+2*Uw2*Aw2*(Tw2-Tz)+q)/Cz
        D(Tw1) ~ (Uw1*Aw1*(Tz-Tw1)+Uw1*Aw1*(To-Tw1))/Cw1
        D(Tw2) ~ (Uw2*Aw2*(Tz-Tw2)+Uw1*Aw1*(To-Tw2))/Cw2
        D(Tr) ~ (Ur*Ar*(Tz-Tr)+Ur*Ar*(To-Tr))/Cr
        D(Wz) ~ (Fsa*(Ws-Wz)+(P/ρa))/Vz]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,100.0)


ev_times = collect(0.0:1.0:100)
condition(u,t,integrator) = t ∈ ev_times
#affect!(integrator) = integrator.u[1] += 5*rand(); print(integrator.p[15])

function affect!(integrator)
    if integrator.p[8] < 0
        integrator.p[8] += (1500*rand())
        integrator.u[1] += (-0.3+0.5*rand())
    else
        integrator.p[8] += (-1500+3000*rand())
        integrator.u[1] += (-0.10+0.08*rand())
    end
    println(integrator.u[1])
end

cb = DiscreteCallback(condition,affect!)

prob = ODEProblem(simpsys,[],tspan,callback=cb,tstops=ev_times)

sol = solve(prob)

f1 = function (p)
    prob1 = remake(prob;p=p)
    sol = solve(prob1,Tsit5();saveat=ev_times)
    return [mean(sol[1,:]), maximum(sol[2,:])]
end

bounds = [[30,100],[0.1,1.0],[1,3],[0.5,2.0],[10,25],[1,3],[1,3],[1,3],[9,9],[12,12],[9,9],[0,1500],[8,40],[70,70],[60,60],[80,80],[36,36],[0,1],[0,2]]

morris_sens = gsa(f1, Morris(), bounds, samples = 100)

fig = Figure(resolution = (600, 400))
fieldnames(typeof(morris_sens))

morris_sens.means_star

Plots.plot(sol)