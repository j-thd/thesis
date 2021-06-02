'''
    Optimize a channel to get the best 

'''

from basic.chamber import hydraulic_diameter_rectangular, radiation_loss, wetted_perimeter_rectangular, velocity_from_mass_flow
import numpy as np
import matplotlib.pyplot as plt

from thermo.prop import FluidProperties
import models.zero_D as zD
import thermo.convection
import basic.IRT as IRT

from thrusters import thruster_data

td = thruster_data.td_Huib_TTH_4_1 # Dictionary holding thruster data for validation and reference


## Design parameters
fp = FluidProperties("water") # Object to get properties of propellant
F = 0.67e-3 # [N] Desired thrust
p_inlet = 4.8e5 # [Pa] Inlet pressure, no pressure drop assumed over chamber
T_inlet = 300 # [K] Inlet temperature
p_back = 0 # [Pa] Into vacuum
h_channel = 100e-6 # [m] Channel depths
#w_channel = 3e-3 # [m] Channel width
AR_exit = 11 # [-] Exit area ratio
T_chamber = 430 # [K]
w_channel_margin = 100e-6 # [m] Padding on each side of the single channel for structural reasons


Nu_func = thermo.convection.Nu_DB # Function used to calculate heat transfer

## Optimization bounds
# Wall temperature (temperature difference with T_chamber)
T_wall_superheat_min = 0 # [K]
T_wall_superheat_max = 400 # [K]
# Channel width
w_channel_min = 1*0.212e-3 # [m] Minimum channel width
w_channel_max = 10.0*0.212e-3 # [m] Maximum channel width

# Saturation temperature
#T_sat = fp.get_saturation_temperature(p=p_inlet)
#print("Saturation temperature: {:3.2f} K".format(T_sat))

## Reference temperature liquid phase
T_bulk_liquid = (T_inlet+T_sat)/2 # [K]
T_bulk_gas = (T_sat+T_chamber)/2 # [K]

## Get required mass flow from IRT engine peformance
ep = zD.engine_performance_from_F_and_T(F_desired=F, p_chamber=p_inlet, T_chamber=T_chamber, AR_exit=AR_exit, p_back=p_back, fp=fp)
m_dot = ep['m_dot'] # [kg/s] Mass flow according to IRT

x_result = zD.minimize_total_power_single_channel(\
    T_wall_superheat_min=T_wall_superheat_min, T_wall_superheat_max=T_wall_superheat_max,\
        w_channel_min=w_channel_min, w_channel_max=w_channel_max, Nu_func=Nu_func,\
            T_inlet=T_inlet, T_chamber=T_chamber, T_ref=T_bulk, p_ref=p_inlet,\
                m_dot=m_dot, h_channel=h_channel, w_channel_margin=w_channel_margin,\
                    emmisivity=emmissivity, fp=fp)

print(x_result)

# Found optimum
T_superheat = x_result.x[0] # [K] 
w_channel = x_result.x[1] # [m]

# Chamber performance following from optimum
T_wall = T_chamber + T_superheat # [K] Wall temperature
A_channel = w_channel * h_channel # [m^2] Channel area
mass_flux = m_dot / A_channel # [kg/m^2] Used as reference for validity of many relations in micro-channel studies
print("A_channel: {}".format(A_channel))
wetted_perimeter = wetted_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Distance of channel cross-section in contact with fluid
print("wetted_perimeter: {}".format(wetted_perimeter))
# Reference length is hydraulic diameter
D_hydraulic = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] Hydraulic diameter

cp = zD.chamber_performance_from_Nu(Nu_func=Nu_func, T_inlet=T_inlet, T_chamber=T_chamber,\
     T_ref=T_bulk, T_wall=T_wall, p_ref=p_inlet, m_dot=m_dot, A_channel=A_channel, L_ref=D_hydraulic, fp=fp)

L_channel = cp['A_heater'] / w_channel # [m]
print("L_channel: {:3.4f} mm".format(L_channel*1e3))
print("mass_flux: {}".format(mass_flux))
print("Mass flow: {:3.4f} mg/s".format(m_dot*1e6))
## Print all results
print(ep)
print(cp)

## Print information for frictional pressure drop guesses
print(" --- FRICTIONAL PRESSURE DROP ESTIMATION ---")
f = 0.05 # [-] Friction factor (emperical based, 0.05 is just a guess)
rho_outlet = fp.get_density(T=T_chamber, p=p_inlet) # [kg/m^3]
u_outlet = velocity_from_mass_flow(A=A_channel, m_dot=m_dot, rho=rho_outlet) # [m/s]
delta_p = f * (L_channel/D_hydraulic) * 0.5 * rho_outlet * u_outlet**2 # [Pa]
print("L/D: {} ".format(L_channel/D_hydraulic))
print("rho_outlet: {}".format(rho_outlet))
print("u_outlet: {}".format(u_outlet))
print("delta_p: {} bar".format(delta_p*1e-5))
print("Relative pressure drop: {}".format(delta_p/p_inlet))