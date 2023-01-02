# HVAC Physical Models

## Zone Model

**definition of a HVAC zone** - An HVAC zone can be defined as a cluster of adjacent offices and/ or spaces often covered by a common air-handling unit (AHU) or air terminal device.

list of analytical models assumptions.


+ The air in the zone is fully mixed, and zone temperature distribution is uniform.
+ The effect of North and South wall on the zone temperature is same. Also, the effect of East and West wall on the zone temperature is same.
+ The ground/floor has no effect on the zone temperature.
+ The density of the air is constant and is not influenced by changing the temperature and humidity ratio of the zone.
+ No pressure losses across the zone and in the mixing section.

### Variables

+ Cz  - Overall thermal capacitance of the zone - 47.1 kJ/C
+ Tz  - Zone Temperature - C
+ Fsa - volume flow rate of supply fan - 0.192 m3/s
+ ρa  - density of air - 1.25 kg/m3 
+ Cpa - specific heat of supply air - 1.005 kJ/kgC
+ Tsa - temperature of supply air - C
+ Uw1 - overall heat transfer coefficient (East and west walls) - 2 w/m2C
+ Uw2 - overall heat transfer coefficient (North and South walls) - 2 w/m2C
+ Ur  - overall heat transfer coefficient (roof) - 1 w/m2C
+ Aw1 - Area (East and West walls) - m2
+ Aw2 - Area (North and South walls) - m2 
+ Ar  - Area (Roof) - m2
+ Tw1 - temperature (East and West walls) - C
+ Tw2 - temperature (North and South walls) - C
+ Tr  - temperature (Roof) - C
+ q   - heat gain for occupants and lights - W

### Model

$$
    C_z\frac{dTz}{dt} = f_{sa}\rho_{a}C_{pa}(T_{sa}-T_z)+2U_{w1}A_{w1}(T_w-Tz)+U_rA_r(T_r-T_z)+2U_{w2}A_{w2}(T_{w2}-T_z)
$$