using  DifferentialEquations, ModelingToolkit, Plots

# variables

# Qc   - Heat transfer rate evaporator
# mᵣ   - Mass transfer rate refrigerent
# h2   - Entalphy (compressor inlet / condensor inlet)
# h6   - Entalphy ( condensor outlet / expansion valve inlet)
# αco  - heat transfer coefficiente (condensor outside)
# αci  - heat transfer coefficiente (condensor inside)
# Aco  - Area (condensor inside)
# Aci  - Area (condensor inside)
# Twc  - Temperature (condensor water)
# Tcw  - Temperature (condensor wall)
# Tc   - Temperature (condensor)

@variables t Twc(t)=0

@parameters αco=0 αci=0 Aco=0 Aci=0 Tcw=0 Tc=0 M=0 C=0

D = Differential(t)

eqs = [
    D(Twc) ~ (αco*Aco*(Tc-Twc)-αci*Aci*(Twc*Tcw))/(M*C)
]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,100.0)


ev_times = collect(0.0:1.0:100)

prob = ODEProblem(simpsys,[],tspan)

sol = solve(prob)

Plots.plot(sol)