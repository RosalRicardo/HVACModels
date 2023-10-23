using  DifferentialEquations, ModelingToolkit, Plots

@variables t Tco(t)=15 Wco(t)=40  

@parameters Cah(t)=4.5 fsw(t)=8.02e-5 ρw(t)=998 Cpw(t)=4.1868 Twi(t)=7 Two(t)=12 UAa(t)=0.04 To(t)=10 fsa(t)=0.192 ρa(t)=1.25 Cpa(t)=1.005 Tm(t)=22 Wm(t)=0.45 Vah(t)=2.88
