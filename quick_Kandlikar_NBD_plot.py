## File to show Nusselt number predicted by Kandlikar

import numpy as np
import matplotlib.pyplot as plt

from thermo.two_phase import Nu_Kandlikar_NBD_dryout
from thermo.prop import FluidProperties
from thermo.convection import Nu_laminar_developed_constant_wall_temp_square, Nu_DB
from basic.chamber import Reynolds_from_mass_flow

fp = FluidProperties("water")

Nu_func_le_laminar = Nu_laminar_developed_constant_wall_temp_square
Nu_func_le_turbulent = Nu_DB
Nu_func_tp = Nu_Kandlikar_NBD_dryout

# Design parameters 
p=5e5 # Pressure through channel [Pa]
m_dot = 25e-6 # [kg/s]
A_channel = 1e-6 # [m^2] Channel cross sections
D_h = 1e-6*np.array((100,200,500)) # [m] Hydraulic diameters

# Mass flux calculated for reference
G = m_dot/A_channel # [kg/(m^2*s)] Mass flux
print("Mass flux: {:3.1f} kg/(m^2*s)".format(G))

# Saturatoin conditions
rho_l = fp.get_liquid_density_at_psat(p_sat=p) # [kg/m^3] Liquid density at saturation point
rho_g = fp.get_vapour_density_at_psat(p_sat=p) # [kg/m^3] Vapour density at saturation point
mu_l = fp.get_liquid_saturation_viscosity(p_sat=p) # [Pa*s] Viscosity at liquid saturatoin point (for Re_le)
Bo = fp.get_Bond_number(p_sat=p, L_ref=D_h) # [-] Bond number








#  Calculate two-phase Nusselt number
x = np.linspace(start=0, stop=1)
D_iter = np.nditer(D_h, flags=['c_index'])
Nu_tp_laminar = np.zeros((np.size(D_h),np.size(x))) # [-] Storing two-phase Nusselt numbers for assumed laminar flow
Nu_tp_turbulent = np.zeros_like(Nu_tp_laminar) # [-] Storing two-phase Nusselt numbers for assumed turbulent flow

for D in D_iter:
    i = D_iter.index
    # Flow parameters for entire flow as a liquid
    Re_le = Reynolds_from_mass_flow(m_dot=m_dot, A_channel=A_channel, L_ref = D, mu=mu_l) # [-] Reynolds as if entire flow is liquid
    print("Re: {}".format(Re_le))
    Pr_le = fp.get_saturation_Prandtl_liquid(p_sat=p) # [-] Prandtl as if entire flow is liquid
    args_le={
        'Re': Re_le,
        'Pr': Pr_le}
    # Calculate Nusselt for entire flow as liquid 
    Nu_le_turbulent = Nu_func_le_turbulent(args=args_le) # [-]
    Nu_le_laminar = Nu_func_le_laminar(args=args_le) # [-]
    # Calculate effect on two-phase Nusselt
    args_tp_laminar = {
        'rho_l':rho_l,
        'rho_g':rho_g,
        'Bo': Bo[i],
        'Nu_le':Nu_le_laminar,
        'Nu_dryout': 0,
        'x': x}
    args_tp_turbulent = {
        'rho_l':rho_l,
        'rho_g':rho_g,
        'Bo': Bo[i],
        'Nu_le':Nu_le_turbulent,
        'Nu_dryout': 0,
        'x': x}

    Nu_tp_laminar[i] = Nu_func_tp(args=args_tp_laminar)
    Nu_tp_turbulent[i] = Nu_func_tp(args=args_tp_turbulent)

fig,axs = plt.subplots(1,2)
axs[0].set_title("Laminar")
axs[1].set_title("Turbulent")

D_iter.reset()
for D in D_iter:
    i = D_iter.index
    axs[0].plot(x,Nu_tp_laminar[i], label="{:3.0f} $\\mu$m".format(D*1e6))
    axs[1].plot(x,Nu_tp_turbulent[i], label="{:3.0f} $\\mu$m".format(D*1e6))

plt.legend(title='Hydraulic diameter $D_h$:')
axs[0].set_ylabel('Nusselt number $Nu_{tp}$ [-]')
axs[0].set_xlabel('Vapour quality $x$ [-]')
axs[0].grid()
axs[1].set_xlabel('Vapour quality $x$ [-]')
axs[1].grid()
plt.suptitle("$Nu$ for Nucleate Boiling $Re<100$ ($p={:1.0f}$ bar, $G={:3.0f}$ kg$\\cdot$m$^{{-2}}$s$^{{-1}}$)".format(p*1e-5,G))
plt.tight_layout(pad=0.5)

plt.show()