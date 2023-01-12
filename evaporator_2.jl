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

@variables t Twe(t)=0

@parameters αeo=0 αei=0 Aeo=0 Aei=0 Tchw=0 Te=0 M=0 C=0

D = Differential(t)

eqs = [
    D(Twe) ~ (αeo*Aeo*(Twe-Te)-αei*Aei*(Tchw*Twe))/(M*C)
]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,100.0)


ev_times = collect(0.0:1.0:100)

prob = ODEProblem(simpsys,[],tspan)

sol = solve(prob)

Plots.plot(sol)