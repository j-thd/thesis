from thermo.prop import FluidProperties
from basic.chamber import ideal_enthalpy_change, required_power

fp = FluidProperties("water")

T_in = 300 # [K]
T_out = 450 # [K]
p_in = 5e5  # [Pa]

m_dot = 1.6667e-6 # [kg/s]

delta_h = ideal_enthalpy_change(T_inlet=T_in, p_inlet=p_in, T_outlet=T_out, p_outlet=p_in, fp=fp) # [J/kg]
Q_dot = required_power(m_dot=m_dot, delta_h=delta_h)

rho_out = fp.get_density(T=T_out, p=p_in)

print("Density out: {:3.3f} kg/m^3".format(rho_out))
print("Required power: {}".format(Q_dot))
