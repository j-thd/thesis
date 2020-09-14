# Script to quickly get some fluid properties at specified temperature and pressure
from thermo.prop import FluidProperties

fp = FluidProperties("water")
T = 886 # [K] Temperature
p = 2.5e5 # [P] Pressure

Pr = fp.get_Prandtl(T=T,p=p) # [-] Prandtl number
mu = fp.get_viscosity(T=T,p=p) # [Pa*s] Viscosity

print("Prandtl number: {:3.2f}".format(Pr))
print("Viscosity {:3.2f} micro*Pa*s".format(mu*1e6))