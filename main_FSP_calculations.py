"""File to calculate the flow similarity parameters from the design parameters established in thrusters.xlsx.
Useful for determining wich relations are interesting and which are not
"""

from basic.chamber import ideal_enthalpy_change, velocity_from_mass_flow, wetter_perimeter_rectangular
from basic.IRT_corrections import hydraulic_diameter
from thermo.prop import FluidProperties

fp = FluidProperties("HEOS::Water")

# Assumed inputs
T_inlet = 300 # [K]
T_outlet = 450 # [K]
p_inlet = 5e5 # [Pa]
h_channel = 100e-6 # [m] Channel height/depth
w_channel = 100e-6 # [m] Channel width
m_dot = 3.6e-6 # [kg/s]

# Hydraulic diameter and area
A_channel = h_channel*w_channel # [m^2] Channel area through which fluid flows
wetted_perimeter  = wetter_perimeter_rectangular(w_channel=w_channel,h_channel=h_channel) # [m] Wetted perimeter
D_h = hydraulic_diameter(A=A_channel, wetted_perimeter=wetted_perimeter) # [m]

print("\n Hydraulic diameter: {:4.4f} um".format(D_h*1e6))

# Determine remaining variables at inlet
rho_inlet = fp.get_density(T=T_inlet,p=p_inlet) # [kg/m^3]
cp_inlet = fp.get_cp(T=T_inlet,p=p_inlet) # [J/(kg*K)] 
u_inlet = velocity_from_mass_flow(m_dot=m_dot, rho=rho_inlet, A=A_channel) # [m/s]
Re_Dh_inlet = fp.get_Reynolds_from_velocity(T=T_inlet, p=p_inlet, L_ref=D_h, u=u_inlet) # [-] Reynolds number at inlet, based on hydraulic diameter
Pr_inlet = fp.get_Prandtl(T=T_inlet,p=p_inlet) # [-] Prandtl number at inlet


print("\n-----     INLET     -----")
print("Density: {:4.8f} kg/m^3".format(rho_inlet))
print("C_p: {:4.8f} kJ/(kg*K)".format(cp_inlet*1e-3))
print("Flow velocity: {:3.4f} m/s".format(u_inlet))
print("Re_Dh: {:4.4f}".format(Re_Dh_inlet))
print("Pr {:4.4f}".format(Pr_inlet))

rho_outlet = fp.get_density(T=T_outlet,p=p_inlet) # [kg/m^3]
cp_outlet = fp.get_cp(T=T_outlet,p=p_inlet) # [J/(kg*K)] 
u_outlet = velocity_from_mass_flow(m_dot=m_dot, rho=rho_outlet, A=A_channel) # [m/s]
Re_Dh_outlet = fp.get_Reynolds_from_velocity(T=T_outlet,p=p_inlet, L_ref=D_h, u=u_outlet) # [-] Re number at outlet, based on hydraulic diameter
Pr_outlet = fp.get_Prandtl(T=T_outlet,p=p_inlet) # [-] Prandtl number at outlet

print("-----     OUTLET    -----")
print("Density: {:4.8f} kg/m^3".format(rho_outlet))
print("C_p: {:4.8f} kJ/(kg*K)".format(cp_outlet*1e-3))
print("Flow velocity: {:3.4f} m/s".format(u_outlet))
print("Re_Dh: {:4.4f}".format(Re_Dh_outlet))
print("Pr {:4.4f}".format(Pr_outlet))

delta_h = ideal_enthalpy_change(T_inlet=T_inlet,p_inlet=p_inlet, T_outlet=T_outlet,p_outlet=p_inlet, fp=fp) # [J/kg] Ideal enthalpy change of propellant
P_delta_h = m_dot * delta_h # [W] Power required to change specific enthalpy for mass flow m_dot

print("-----     TOTAL    -----")
print("P_delta_h: {:3.4f} W".format(P_delta_h))
