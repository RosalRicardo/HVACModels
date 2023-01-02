using DifferentialEquations, ModelingToolkit, Plots

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

@variables t Tz(t)=20

@parameters Cz=47.1 Fsa=0.192  ρa=1.25 Cpa=1.005 Tsa=15 Uw1=2 Uw2=2 Ur=1 Aw1=9 Aw2=12 Ar=9 Tw1=18 Tw2=17 Tr=25 q=300

D = Differential(t)

eqs = [D(Tz) ~ (Fsa*ρa*Cpa*(Tsa-Tz)+2*Uw1*Aw1*(Tw1-Tz)+Ur*Ar*(Tr-Tz)+2*Uw2*Aw2*(Tw2-Tz)+q)/Cz]

@named sys = ODESystem(eqs,t)

simpsys = structural_simplify(sys)

tspan = (0.0,10.0)
prob = ODEProblem(simpsys,[],tspan)

sol = solve(prob)

plot(sol)