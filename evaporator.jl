using  DifferentialEquations, ModelingToolkit, Plots

# variables

# Qₑ   - Heat transfer rate evaporator
# mᵣ   - Mass transfer rate refrigerent
# h1   - Entalphy (evaporator outlet / compressor inlet)
# h6   - Entalphy (expansion valve exit / evaporator inlet)
# αeo  - heat transfer coefficiente (evaporator outside)
# αei  - heat transfer coefficiente (evaporator inside)
# Aeo  - Area (evaporator inside)
# Aei  - Area (evaporator inside)
# Tchw - Temperature (evaporator coolant water)
# Twe  - Temperature (evaporator wall)
# Twe  - Temperature (evaporator)

@variables t Qₑ(t)=0 mᵣ(t)=0

@parameters h1=0 h6=0 αeo=0 αei=0 Aeo=0 Aei=0 Tchw=0 Twe=0 Te=0

D = Differential(t)

eqs = [
    D(Qₑ) ~ αei*Aei*(Tchw-Twe)
    D(mᵣ) ~ (αeo*Aeo*(Tchw-Twe))
]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,100.0)


ev_times = collect(0.0:1.0:100)

prob = ODEProblem(simpsys,[],tspan)

sol = solve(prob)

Plots.plot(sol)