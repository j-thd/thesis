## File to plot Nusselt number relations for the report

import numpy as np
import matplotlib.pyplot as plt

from thermo.prop import FluidProperties
import thermo.convection as conv
import basic.chamber

# Object to retrievere propellant properties from
fp = FluidProperties("water")

# Reference state for calculations
p = 5e5 # [Pa]
T_liquid = 300 # [K]
T_gas = 500 # [K]
m_dot = 1e-6 # [kg/s] For Kandlikar relation.
#h_c = 100e-6 # [m] Channel depth/height

# Density is needed to reverse calculate the reference length for Re, so it can be applied in other flow similarity parameters
rho_liquid = fp.get_density(T=T_liquid, p=p)
mu_liquid = fp.get_viscosity(T=T_liquid, p=p)

# Range of Reynolds numbers to consider
Reynolds = np.logspace(start=0, stop=2.8)
D_hydraulic = np.logspace(start=-4.93, stop=-2.93)
A_channel = D_hydraulic**2
Reynolds_2 = fp.get_Reynolds_from_mass_flow(T=T_liquid,p=p,L_ref=D_hydraulic, m_dot=m_dot, A=A_channel)
# Range of hydraulic diameter that fit this Reynolds number
#w_channel = 2*m_dot/( Reynolds_100 * mu_liquid ) - h_c # [m]
#A_channel = h_c * w_channel # [m^2]
#D_hydraulic = basic.chamber.hydraulic_diameter_rectangular(w_channel=w_channel,h_channel=h_c) # [m]

#Reynolds_check = fp.get_Reynolds_from_mass_flow(T=T_liquid, p=p, L_ref=D_hydraulic, m_dot=m_dot, A=A_channel) # [-]
Bond = fp.get_Bond_number(p_sat=p, L_ref=D_hydraulic) # [m]
# Prandtl number in both reference states
Pr_liquid = fp.get_Prandtl(T=T_liquid,p=p) # [-]
Pr_gas = fp.get_Prandtl(T=T_gas,p=p) # [-]
# Conductivity in reference state
conductivity_liquid = fp.get_thermal_conductivity(T=T_liquid, p=p)




Nu_DB_liquid = np.zeros_like(Reynolds)
Nu_DB_gas = np.zeros_like(Reynolds)  
Nu_Kandlikar_liquid = np.zeros_like(Reynolds)
#Nu_Kandlikar_gas = np.zeros_like(Reynolds)
# Calculate Nusselt number for different relations

# Iterare over Reynolds numbers
# print(Reynolds)
# iter_Re = np.nditer(Reynolds, flags=['c_index'])
# for Re in Reynolds:

args_liquid = {'Re': Reynolds, 'Pr':Pr_liquid}
args_gas = {'Re': Reynolds, 'Pr':Pr_gas}
args_Kandlikar = {'Bo': Bond, 'D_hydraulic': D_hydraulic, 'conductivity': conductivity_liquid}
Nu_DB_liquid = conv.Nu_DB(args=args_liquid)
Nu_DB_gas = conv.Nu_DB(args=args_gas)
Nu_FDL_wall  = conv.Nu_laminar_developed_constant_wall_temp_square(args={})
Nu_FDL_flux  = conv.Nu_laminiar_developed_constant_heat_flux_square(args={})
Nu_Kandlikar_liquid = conv.Nu_Kandlikar_NDB_Re_low_sat_gas_constant_wall_temp_square_water(args=args_Kandlikar)


plt.semilogx(Reynolds, Nu_DB_liquid, label="DB ($T={:3.0f}$ K, $Pr={:1.2f}$)".format(T_liquid, Pr_liquid))
plt.semilogx(Reynolds, Nu_DB_gas, label="DB ($T={:3.0f}$ K, $Pr={:1.2f}$)".format(T_gas, Pr_gas))

plt.hlines(Nu_FDL_wall, xmin=Reynolds[0], xmax=Reynolds[-1], colors='r', linestyle='dashed', label="FDL square ($T_w = const.$)")
plt.hlines(Nu_FDL_flux, xmin=Reynolds[0], xmax=Reynolds[-1], colors='r', linestyle='dotted', label="FDL square ($q = const.$)")


plt.xlabel('Reynolds number $Re$ [-]')
plt.ylabel('Nusselt number $Nu$ [-]')
plt.title('Nusselt number relations for square channel ($p={:1.0f}$ bar, $\dot{{m}}={:1.0f}$ mg/s )'.format(p*1e-5, m_dot*1e6))
plt.grid()
plt.legend()

plt.figure()
plt.semilogx(Reynolds_2, Nu_Kandlikar_liquid, label="Kandlikar Re < 100")
plt.xlabel('Reynolds number $Re$ [-]')
plt.ylabel('Nusselt number $Nu$ [-]')
plt.title('Nusselt number relations for square channel ($p={:1.0f}$ bar, $\dot{{m}}={:1.0f}$ mg/s )'.format(p*1e-5, m_dot*1e6))
plt.legend()
plt.grid()

plt.show()