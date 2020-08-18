import math
from models.zero_D import DB_IRT, DB_Rajeev
import thermo.prop

# Chamber input values
p_inlet = 5e5 # [Pa] Inlet pressure (preferably 5 bar)
T_inlet = 295.86 # [K] Inlet temperature (usually room temperature)
T_wall = 1000 # [K] Wall temperature
L_channel = 100e-3 # [m] Channel length
w_channel = 1e-3 # [m] Channel width
h_channel = 1e-3 # [m] Channel height (also called channel depth)
u = 1 # [m/s] Flow velocity


# Nozzle design values
w_throat = 1.74e-5 # [m] Throat width
h_throat = 8.1e-5 # [m] Throat height (sometimes known as channel depth)
throat_roc = 1e-6 # [m] Throat radius of curvature

# Exit parameters
AR_exit = 30 # [m] Area ratio (exit/throat)
p_back = 0 # [Pa] Back pressure
divergence_half_angle = math.radians(20.5) # [rad] Divergence angle of nozzle cone

# Fluid of choice
fp = thermo.prop.FluidProperties("nitrogen")

#DB_Rajeev(T_inlet=T_inlet,T_wall=T_wall,p_inlet=p_inlet,u=u,L_channel=L_channel, w_channel=w_channel,h_channel=h_channel,w_throat=w_throat, \
    #h_throat=h_throat,throat_roc=throat_roc,AR_exit=AR_exit,p_back=p_back,divergence_half_angle=divergence_half_angle,fp=fp)

## 

F_desired = 0.46e-3 # [N] Desired thrust  
p_inlet = 1.92e5 # [Pa] Inlet pressure
#h_throat = 8.1e-5 # [m] Throat height
#w_throat = 1.74e-5 # [m] Throat width
AR_exit = 30
p_back = 0
T_chamber = 500 # [K] Chamber temperature
T_inlet = 295.86 # [K] Inlet temperature (usually room temperature)
L_channel = 100e-3 # [m] Channel length

results = DB_IRT(F_desired=F_desired,p_inlet=p_inlet,T_inlet=T_inlet ,T_chamber=T_chamber, L_channel=L_channel, w_channel=w_channel, h_channel=h_channel,\
    w_throat=w_throat, h_throat=h_throat, AR_exit=AR_exit, p_back=p_back, fp=fp)