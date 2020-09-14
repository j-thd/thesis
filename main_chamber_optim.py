# File to optimize channel based on simple Nusselt number
from thermo.convection import Nu_DB, Stanton_from_Nusselt_func_and_velocity
import numpy as np
from thermo.prop import FluidProperties
from basic.chamber import h_conv_from_Stanton, hydraulic_diameter_rectangular, ideal_enthalpy_change, radiation_loss, required_heater_area, required_power, velocity_from_mass_flow
import matplotlib.pyplot as plt


# Given parameters
fp = FluidProperties("water") # Object to get water properties from
T_inlet = 300 # [K] Temperature at inlet (room temperature)
p_inlet = 5e5 # [Pa] Pressure at inlet
T_chamber = 900 # [K] Desired exit temperature
T_wall = 1000 # [K] Given wall temperature
h_channel = 100e-6 # [m] Channel depth
m_dot = 1e-6 # [kg/s] Mass flow
# Used Nusselt relation
Nu_func = Nu_DB


# Which channel width gives the least power usage?
w_channel_min = 1e-3 # [m] Minimum channel width
w_channel_max = 5e-3 # [m] Maximum channel width
w_channel = np.linspace(start=w_channel_min, stop=w_channel_max) # [m] Range of channel widths to evaluate

# Calculate some basic channel geometry
A_channel = w_channel*h_channel # [m^2] Cross-sectional area of channel
Dh_channel = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter of rectangular channels

# Calculate some basic flow parameter at the inlet
rho_inlet = fp.get_density(T=T_inlet, p=p_inlet) # [kg/m^3]
u_inlet = velocity_from_mass_flow(m_dot=m_dot, rho=rho_inlet, A=A_channel) # [m/s] Flow velocity inside channel

# Calculate heat transfer paramaters
T_bulk = (T_inlet + T_chamber) / 2 # [K] Bulk temperature as reference
Stanton = Stanton_from_Nusselt_func_and_velocity(Nu_func=Nu_func, u=u_inlet, T_ref=T_bulk, p_ref=p_inlet, L_ref=Dh_channel, fp=fp) # [-] Stanton number
h_conv = h_conv_from_Stanton(Stanton=Stanton, u=u_inlet, T_ref=T_bulk, p_ref=p_inlet, fp=fp) # [W/(m^2*K)] Convective heat transfer coefficient
delta_h = ideal_enthalpy_change(T_inlet=T_inlet, p_inlet=p_inlet, T_outlet=T_chamber, p_outlet=p_inlet, fp=fp) # [J/kg]
Q_dot = required_power(m_dot=m_dot,delta_h=delta_h) # [W] Required power to increase enthalpy to delta_h
# Resulting heater area
A_heater = required_heater_area(Q_dot=Q_dot, h_conv=h_conv, T_wall=T_wall, T_ref=T_bulk) # [m^2] Required heater area
L_channel = A_heater/w_channel # [m] Final length of channel for the give w_channel
# Radiaton loss
P_radiation_loss = radiation_loss(T=T_wall, A=A_heater, emmisivity=1) # [W] Heat radiation to environment
P_total = Q_dot+P_radiation_loss # [W] Total required power

# Print some results
print("Density inlet: {:4.2f} kg/m^3".format(rho_inlet))
print("Bulk Prandtl {:2.2f}".format(fp.get_Prandtl(T=T_bulk,p=p_inlet)))
print("Q_dot: {}".format(Q_dot))

#plt.plot(w_channel*1e3, fp.get_Reynolds_from_velocity(T=T_bulk,p=p_inlet,L_ref=Dh_channel,u=u_inlet))
plt.plot(w_channel*1e3, A_heater*1e6)
plt.xlabel("Channel width $w_c$ [mm]")
plt.ylabel("Heater area [mm^2]")
plt.grid()
plt.title("heater area")

plt.figure()
plt.plot(w_channel*1e3, L_channel*1e3)
plt.xlabel("Channel width $w_c$ [mm]")
plt.ylabel("Channel length [mm]")
plt.grid()
plt.title("Channel length")

plt.figure()
plt.plot(w_channel*1e3, fp.get_Reynolds_from_velocity(T=T_bulk,p=p_inlet,L_ref=Dh_channel,u=u_inlet))
plt.xlabel("Channel width $w_c$ [mm]")
plt.ylabel("Reynolds [-]")
plt.grid()
plt.title("Reynolds")

plt.figure()
plt.plot(w_channel*1e3, P_total)
plt.xlabel("Channel width $w_c$ [mm]")
plt.ylabel("Power [W]")
plt.grid()
plt.title("Power")



plt.show()