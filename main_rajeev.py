import math
import basic.IRT_corrections
import thermo.prop
import thrusters.thruster_data
from physical_constants import g0

td = thrusters.thruster_data.Cen2010_6 # {dict} with thruster data


# Chamber input values
p_chamber = td['p_inlet'] # [Pa]
T_chamber = td['T_chamber_guess'] # [K]

print("Chamber pressure: {:1.2f} bar".format(p_chamber*1e-5))
print("Chamber temperature: {:3.2f} K".format(T_chamber))

# Nozzle design values
w_throat = td['w_throat'] # [m] Throat width
h_throat = td['h_channel'] # [m] Throat height (sometimes known as channel depth)
throat_roc = 1e-6 # [m] Throat radius of curvature

# Exit parameters
AR_exit = td['AR_exit']# [m] Area ratio (exit/throat)
p_back = td['p_back'] # [Pa] Back pressure
divergence_half_angle = td['divergent_half_angle'] # [rad] Divergence angle of nozzle cone

# Fluid of choice
fp = thermo.prop.FluidProperties(td['propellant'])
is_cold_flow = False # True if the flow is not heated

basic.IRT_corrections.Rajeev_complete(p_chamber=p_chamber,T_chamber=T_chamber,w_throat=w_throat,h_throat=h_throat\
    , throat_roc=throat_roc, AR_exit=AR_exit, p_back=p_back,divergence_half_angle=divergence_half_angle, fp=fp, is_cold_flow=is_cold_flow)

print("\n === EXPERIMENTAL VALUES === ")
print(" Thrust: {:3.2f} mN".format(td['F']*1e3))
print(" Mass flow: {:3.2f} mg/s".format(td['m_dot']*1e6))
Isp_experimental = td['F']/td['m_dot']/g0
print(" Isp: {:3.2f} s".format(Isp_experimental))
