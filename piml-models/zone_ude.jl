using  Flux, CairoMakie, DifferentialEquations, ModelingToolkit, Plots, GlobalSensitivity, Statistics

# variables
# Tz  - Zone Temperature - C
# Tw1 - temperature (East and West walls) - C
# Tw2 - temperature (North and South walls) - C
# Tr  - temperature (Roof) - C

# parameters
# 01 Cz  - Overall thermal capacitance of the zone - 47.1 kJ/C
# 02 Fsa - volume flow rate of supply fan - 0.192 m3/s
# 03 ρa  - density of air - 1.25 kg/m3 
# 04 Cpa - specific heat of supply air - 1.005 kJ/kgC
# 05 Tsa - temperature of supply air - C
# 06 Uw1 - overall heat transfer coefficient (East and west walls) - 2 w/m2C
# 07 Uw2 - overall heat transfer coefficient (North and South walls) - 2 w/m2C
# 08 Ur  - overall heat transfer coefficient (roof) - 1 w/m2C
# 09 Aw1 - Area (East and West walls) - m2
# 10 Aw2 - Area (North and South walls) - m2 
# 11 Ar  - Area (Roof) - m2
# 12 T₀  - Outside Air Temperature - C
# 13 Cw1 - 
# 14 Cw2 - 
# 15 Cr  - 
# 16 Vz  - 
# 17 Ws  - 
# 18 P   - 

# 12 q   - heat gain for occupants and lights - W

u0 = Float32[40.0;20.0;20.0;25.0;0.5;3000]
tspan = (0.0f0,100.0f0)

ann = Chain(Dense(6,10,tanh), Dense(10,1))
p1,re = Flux.destructure(ann)
ps = Flux.params(p3)

# @parameters Cz=47.1e3 Fsa=0.192  ρa=1.25 Cpa=1.005 Tsa=8 Uw1=2 Uw2=2 Ur=1 Aw1=9 Aw2=12 Ar=9 q=3000 To=21 Cw1=70 Cw2=60 Cr=80 Vz=36 Ws=0.02744 P=0.08 weights=p1 ann=re

p2 = Float32[47000.0;0.192;1.25;1.005;8.0;2.0;2.0;1.0;9.0;12.0;9.0;21.0;70.0;60.0;80.0;36.0;0.02744;0.08]
p3 = [p1;p2]
ps = Flux.params(p3)

function dudt_(du,u,p,t)
    Tz, Tw1, Tw2, Tr, Wz, load = u
    du[1] = (p[82]*p[84]*p[85]*(p[86]-Tz)+2*p[87]*p[90]*(Tw1-Tz)+p[89]*p[92]*(Tr-Tz)+2*p[88]*p[91]*(Tw2-Tz)+load)/p[82]
    du[2] = (p[87]*p[90]*(Tz-Tw1)+p[87]*p[90]*(p[93]-Tw1))/p[94]
    du[3] = (p[88]*p[91]*(Tz-Tw2)+p[88]*p[91]*(p[93]-Tw2))/p[95]
    du[4] = (p[89]*p[92]*(Tz-Tr)+p[89]*p[82]*(p[93]-Tr))/p[96]
    du[5] = (p[83]*(p[88]-Wz)+(p[99]/p[84]))/p[97]
    du[6] = re(p[1:81])(u)[1]
end

prob = ODEProblem(dudt_,u0,tspan,p3)
concrete_solve(prob,Tsit5(),u0,p3,abstol=1e-8,reltol=1e-6)

function predict_adjoint()
    Array(concrete_solve(prob,Tsit5(),u0,p3,saveat=0.0:0.1:100.0))
  end
  loss_adjoint() = sum(abs2,x-1 for x in predict_adjoint())
  loss_adjoint()
  
  data = Iterators.repeated((), 300)
  opt = ADAM(0.01)
  iter = 0
  cb = function ()
    global iter += 1
    if iter % 50 == 0
      display(loss_adjoint())
      display(plot(solve(remake(prob,p=p3,u0=u0),Tsit5(),saveat=0.1),ylim=(0,6)))
    end
  end
  
  # Display the ODE with the current parameter values.
  cb()
  
  Flux.train!(loss_adjoint, ps, data, opt, cb = cb)

