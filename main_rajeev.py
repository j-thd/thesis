import math
import basic.IRT_corrections
import thermo.prop

# Chamber input values
p_chamber = 1.92e5 # [Pa]
T_chamber = 295.86 # [K]

# Nozzle design values
w_throat = 1.74e-5 # [m] Throat width
h_throat = 8.1e-5 # [m] Throat height (sometimes known as channel depth)
throat_roc = 1e-6 # [m] Throat radius of curvature

# Exit parameters
AR_exit = 30 # [m] Area ratio (exit/throat)
p_back = 0 # [Pa] Back pressure
divergence_angle = math.radians(20.5) # [rad] Divergence angle of nozzle cone

# Fluid of choice
fp = thermo.prop.FluidProperties("nitrogen")
is_cold_flow = True # True if the flow is not heated

basic.IRT_corrections.Rajeev_complete(p_chamber=p_chamber,T_chamber=T_chamber,w_throat=w_throat,h_throat=h_throat\
    , throat_roc=throat_roc, AR_exit=AR_exit, p_back=p_back,divergence_half_angle=divergence_angle, fp=fp, is_cold_flow=is_cold_flow)