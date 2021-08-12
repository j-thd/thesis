## File to contain the 1D- model calculations. Heat flux relations, etc... come from elsewhere

from thermo.prop import FluidProperties
from basic.chamber import inlet_manifold_width, velocity_from_mass_flow, Reynolds_from_mass_flow, required_power, required_heater_area, Froude_number, delta_enthalpy_per_section,\
    wetted_perimeter_rectangular, hydraulic_diameter_rectangular
from thermo.convection import heat_transfer_coefficient_from_Nu
import thermo.two_phase as tp
import numpy as np

## Two functions for LIQUID PHASE 

def prepare_single_phase_liquid(T_inlet, steps, p_ref, m_dot, fp: FluidProperties):
    """ Prepare numpy arrays for calculating channel length in a liquid single-phase section of a channel.
    NOTE: This is done to avoid recalculating arrays that are not dependent on channel geometry, therefore speeding up optimizations.\
        After all, during optimization the geometry is what varies.\
        Also it also ensure that temperature endpoint and enthalpy cleanly match with saturation temperature in the correct phase

    Args:
        T_inlet (K): Inlet temperature
        steps (-): Amount steps of dT taken to reach saturation temperature T_sat (dT = (T_sat-T_inlet)/2)
        p_ref (Pa): Pressure assumed constant along channel, equal to inlet pressure
        m_dot (kg/s): Mass flow
        fp (FluidProperties): Object to access propellant properties with
    """
    T_sat = fp.get_saturation_temperature(p=p_ref) # [K] Saturation temperature
    assert ( T_inlet < T_sat) # Check input
    assert (steps > 1)

    # Temperature and other intermediate variable in channel section i=0...n
    T, dT = np.linspace(start=T_inlet, stop=T_sat, num=steps,retstep=True) # [K] Temperature T_i (also returns steps between sections)
    # The reference temperature for heat transfer calculations
    # The first value [0] should not be important. The heat transfer calculated at i is between i-1 and i
    #  So, from T[i-1] to T[i]. So, if there reference temperature is the average dT/2 must SUBTRACTED
    #T_ref = T - dT/2 # [K] Reference temperature for heat transfer calculations

    ## Get all thermodynamic values that can be precalculated
    # NOTE: all last values must be replaced with the correct values for the saturated liquid state
    # Before the values are replaced, sometimes an error is thrown because the values are close to the saturation point
    # That, or NaNs and infinites show up. This shouldn't be a problem, unless the second-to-last points also start getting close to the saturation point
    
    # Enthalpy 
    h = fp.get_enthalpy(T=T, p=p_ref) # [J/kg] Enthalpy
    h[-1] = fp.get_saturation_enthalpy_liquid(p=p_ref) # [J/kg] Saturation enthalpy at T_n = T_sat
    # Heating power required in section to increase temp by dT. Use enthalpy difference
    delta_h = delta_enthalpy_per_section(h=h) # [J/kg] Enthalpy difference per section
    Q_dot = required_power(m_dot=m_dot, delta_h=delta_h) # [W]

    # Density
    rho = fp.get_density(T=T, p=p_ref) # [kg/m^3] Density
    rho[-1] = fp.get_liquid_density_at_psat(p_sat=p_ref) # [kg/m^3] Saturation density
    # Prandtl number
    Pr = fp.get_Prandtl(T=T, p=p_ref) # [-] Prandtl number
    Pr[-1] = fp.get_saturation_Prandtl_liquid(p_sat=p_ref) # [-] Saturation Prandtl
    # Thermal conductivity 
    kappa = fp.get_thermal_conductivity(T=T, p=p_ref) # [W/(m*K)] Conductivity
    kappa[-1] = fp.get_liquid_saturation_conductivity(p_sat=p_ref) # [W/(m*K)] Saturation conductivity
    # Viscosity
    mu = fp.get_viscosity(T=T, p=p_ref) # [Pa*s] Viscosity
    mu[-1] = fp.get_liquid_saturation_viscosity(p_sat=p_ref) # [Pa*s] Saturation viscosity
    return {\
        "T":T, # [K]
        "dT": dT, # [K]
        "rho": rho, # [kg/m^3]
        "h": h, # [J/kg]
        "Q_dot": Q_dot, # [W]
        "Pr": Pr, # [-]
        "kappa": kappa, # [W/(m*K)]
        "mu": mu, # [Pa*s]
        }

def calc_channel_single_phase(T, Q_dot, rho, Pr, kappa, mu, p_ref, m_dot, T_wall, D_hydr, wetted_perimeter, A_channel, Nu_func, fp: FluidProperties):
    """ Calculate channel length and intermediate variables by using a constant temperature step DT
        Using DT should give the best accuracy for computation time
        Must be a single-phase channel.

        NOTE: for optimization, this calculation should be preceded by a prepatory calculations, which really only needed to be once (enthalpy, etc..)

    Args:
        T_0 (K): Temperature  at channel entrance
        T_n (K): Temperature at channel exit 
        steps (-): Temperature steps into which channel is divided (dT = (T_n - T_0)/steps)
        p_ref (Pa): Assume pressure is constant through channel, p_ref = p_inlet
        m_dot (kg/s): Mass flow
        wetted_perimeter (m): Channel perimeter through which heat can be conducted
        A_channel (m^2): Area through which fluid flows
        fp (FluidProperties): Object to access fluid properties with
    """

    u = velocity_from_mass_flow(A=A_channel, m_dot=m_dot, rho=rho) # [m/s]
    # Reynold's is not calculated with FluidProperties to avoid the issues at the saturation point (which sohuld have been handled in the prepare_single_phase function)
    Re = Reynolds_from_mass_flow(m_dot=m_dot, A_channel=A_channel, L_ref=D_hydr, mu=mu) # [-]
    # The Nusselt number is calculated with pre-calculated Prandtl number, and the geometry-specific Reynolds number
    Nu = Nu_func(args=  { 
        "Re": Re,
        "Pr": Pr,
        })
    # Convective heat transfer parameter
    h_conv = heat_transfer_coefficient_from_Nu(Nu=Nu, kappa=kappa, L_ref=D_hydr)
    # Calculate the heat flux in a channel section
    A_heating = required_heater_area(Q_dot=Q_dot, h_conv=h_conv, T_wall=T_wall, T_ref=T) # [m^2] Heating area required to deliver Q_dot
    delta_L = A_heating / wetted_perimeter # [m] Length of channel sections to obtain A_heating
    L = np.cumsum(delta_L)
    
    # Return a dictionary with values of interest
    return {
        # All values are for given temperature T
        'L': L, # Distance along channel at temperature T
        'Nu': Nu, # [-] Nusselt number (based on L_ref)
        'Re': Re, # [-] Reynolds number (based on L_ref)
        'u': u, # [m/s] Flow velocity
        'h_conv': h_conv, # [W/(m*K)] Convective heat transfer parameter 
        'delta_L': delta_L, # [m] Length of each section
    }

## Two functions for multi-phase saturation phase from x=0 to x=1
# Homogenous model in which both phases have same velocity u_g = u_l

def prepare_homogenous_transition(p, m_dot, steps, fp: FluidProperties):
    x = np.linspace(start=0, stop=1, num=steps) # [-] Vapour quality range
    ## NOTE: subscript sat for sat has been dropped for readability
    # Calculate saturation parameters at edges
    T_sat = fp.get_saturation_temperature(p=p) # [K] Saturation temperature
    rho_l = fp.get_liquid_density_at_psat(p_sat=p) # [kg/m^3]
    rho_g = fp.get_vapour_density_at_psat(p_sat=p) # [kg/m^3] Gas saturation density

    # Void fraction is precalculated because it allow for simple evaluation of velocity when geometry changes
    alpha = tp.homogenous_void_fraction(x=x, rho_g=rho_g, rho_l=rho_l) # [-] Void fraction
    rho = tp.mixture_density(alpha=alpha, rho_g=rho_g, rho_l=rho_l) # [kg/m^3] Mixture density of two-phase flow

    # Mean viscosity has no obvious way to be calculated and as such, a relation must simply be chosen [10.42] from Carey2008 is used.
    mu_l = fp.get_liquid_saturation_viscosity(p_sat=p) # [Pa*s]
    mu_g = fp.get_gas_saturation_viscosity(p_sat=p) # [Pa*s]
    mu = tp.mean_viscosity(mu_g=mu_g, mu_l=mu_l, rho_l=rho_l, rho_g=rho_g, x=x) # [Pa*s]

    # Thermal conductivity at saturation
    kappa_l = fp.get_liquid_saturation_conductivity(p_sat=p) # [W/(m*K)]
    kappa_g = fp.get_gas_saturation_conductivity(p_sat=p) # [W/(m*K)]
    # Mean conductivity
    kappa = tp.mean_conductivity(kappa_g=kappa_g, kappa_l=kappa_l, rho_l=rho_l, rho_g=rho_g, x=x) # [W/(m^2*K)]
    # Prandtl numbers at saturation,
    Pr_l = fp.get_saturation_Prandtl_liquid(p_sat=p) # [-]
    Pr_g = fp.get_saturation_Prandtl_gas(p_sat=p) # [-]
    # Mean Prandtl
    Pr = tp.mean_Prandtl(Pr_g=Pr_g, Pr_l=Pr_l, rho_l=rho_l, rho_g=rho_g, x=x) # [-]



       


    # Saturation enthalpies
    h_sat_liquid = fp.get_saturation_enthalpy_liquid(p=p) # [J/kg]
    h_sat_gas = fp.get_saturation_enthalpy_gas(p=p) # [J/kg]
    # Enthalpy as function of vapour quality x
    h = h_sat_liquid + (h_sat_gas-h_sat_liquid) * x # [J/kg] Saturation enthalpy as flow quality increases
    delta_h = delta_enthalpy_per_section(h=h) # [J/kg] Enthalpy difference per section
    Q_dot = required_power(m_dot=m_dot, delta_h=delta_h) # [W] Heating power required to increase enthalpy in each sections

    return {
        'x': x,
        'alpha': alpha,
        'T_sat': T_sat,
        'rho': rho,
        'rho_l': rho_l,
        'rho_g': rho_g,
        'mu': mu,
        'mu_l': mu_l,
        'mu_g': mu_g,
        'Pr_l': Pr_l,
        'Pr_g': Pr_g,
        'Pr': Pr,
        'kappa_l': kappa_l,
        'kappa_g': kappa_g,
        'kappa': kappa,
        'h': h,
        'Q_dot': Q_dot,
    }

def calc_homogenous_transition(p_sat, x, alpha, T_sat, rho_l, rho_g, rho, m_dot, mu_l, mu, Pr, Pr_l, kappa, kappa_l, Q_dot, T_wall, D_hydr, wetted_perimeter, A_channel, Nu_func_tp, Nu_func_le, Nu_func_dryout, fp: FluidProperties):
    
    # Calculate flow velocity 
    # Homogenous assumption means that liquid and gas velocity are equal, so one of two must be calculated.
    # Alpha (void fraction) has already been calculated to satisfy this condition
    # For both velocities there will be some issues with a singularity at (x=0, alpha =0) or (x=1, alpha=1), but the overal density rho does not have this problem
    #u_l = (m_dot / A_channel) * ( 1 - x ) /  (rho_l * ( 1- alpha )) # [m/s]
    #u_g = (m_dot / A_channel) *     x     /  (rho_g *     alpha) # [m/s]
    u = (m_dot / A_channel) / rho # [m/s]
    
    ## Calculate flow similarity/two-phase parameters that can only be known as soon as the geometry is known
    ## NOTE: mu is some form of weighted average and there are several ways to weight it (prepare_ function has the method)
    Re = Reynolds_from_mass_flow(m_dot=m_dot, A_channel=A_channel, L_ref=D_hydr, mu=mu) # [-]

    # Bond number
    Bo = fp.get_Bond_number(p_sat=p_sat,L_ref=D_hydr) # [-] Bond number

    # Nusselt number for entire flow as liquid (often used for two-phase relations)
    # NOTE: _le is often used as subscript, but _l might sometimes mean something else, such as G(1-x) * L / (A * mu_l). _le does NOT scale with (1-x) and is G * L / ( A * mu_l ) 
    Re_le = Reynolds_from_mass_flow(m_dot=m_dot, A_channel=A_channel, L_ref=D_hydr, mu=mu_l) # [-] Use viscosity of liquid phase to get Re_le
    args_le = { 'Re': Re_le,
                'Pr': Pr_l}
    
    Nu_le = Nu_func_le(args=args_le) # [-] Nusselt number as if entire flow was liquid.
    args_dryout = { 
        'Re': Re,
        'Pr': Pr,
        'kappa': kappa,
        'kappa_l': kappa_l,
    }
    Nu_dryout = Nu_func_dryout(args=args_dryout)

    ## Actual Nusselt number
    args = {    'Re': Re,
                'Pr': Pr_l,
                'Bo': Bo,
                'rho_g': rho_g,
                'rho_l': rho_l,
                'x': x,
                'Nu_le': Nu_le,
                'Nu_dryout': Nu_dryout,
                }
    Nu = Nu_func_tp(args=args) # [-] Two-phase Nusselt number

    # The Nusselt number is calculated in relation to liquid thermal conducitvity only
    h_conv = heat_transfer_coefficient_from_Nu(Nu=Nu, kappa=kappa_l, L_ref=D_hydr) # [W/(m^2 * K)] Heat transfer coefficient
    A_heater = required_heater_area(Q_dot=Q_dot, h_conv=h_conv, T_wall=T_wall, T_ref=T_sat) # [m^2] Required heater area to heat up section with Q_dot
    delta_L = A_heater / wetted_perimeter # [m] Required channel section length for heating area
    L = np.cumsum(delta_L) # [m] Cumulative channel length of two-phase section


    return {
        'L': L,
        'delta_L': delta_L,
        'u': u,
        'rho': rho,
        'Re': Re,
        'Nu': Nu,
        'Nu_le': Nu_le,
        'h_conv': h_conv,
    }

## Just one extra function needed for the gas phase, as inputs should be the same as single-phase function for liquid

def prepare_single_phase_gas(T_outlet, steps, p_ref, m_dot, fp: FluidProperties):

    T_sat = fp.get_saturation_temperature(p=p_ref) # [K] Saturation temperature
    assert (T_outlet > T_sat)
    assert (steps > 1)

    # Temperature and other intermediate variable in channel section i=0...n
    T, dT = np.linspace(start=T_sat, stop=T_outlet, num=steps, retstep=True) # [K] Temperature T_i
    # The reference temperature for heat transfer calculations
    # The first value [0] should not be important. The heat transfer calculated at i is between i-1 and i
    #  So, from T[i-1] to T[i]. So, if there reference temperature is the average dT/2 must SUBTRACTED
    #T_ref = T - dT/2 # [K] Reference temperature for heat transfer calculations

    ## Get all thermodynamic values that can be precalculated
    # NOTE: all first values must be replaced with the correct values for the saturated gas state
    # Before the values are replaced, sometimes an error is thrown because the values are close to the saturation point
    # That, or NaNs and infinites show up. This shouldn't be a problem, unless the second-to-last points also start getting close to the saturation point
    
    # Enthalpy 
    h = fp.get_enthalpy(T=T, p=p_ref) # [J/kg] Enthalpy
    h[0] = fp.get_saturation_enthalpy_gas(p=p_ref) # [J/kg] Saturation enthalpy at T_n = T_sat
    # Heating power required in section to increase temp by dT. Use enthalpy difference
    delta_h = delta_enthalpy_per_section(h=h) # [J/kg] Enthalpy difference per section
    Q_dot = required_power(m_dot=m_dot, delta_h=delta_h) # [W]

    # Density
    rho = fp.get_density(T=T, p=p_ref) # [kg/m^3] Density
    rho[0] = fp.get_vapour_density_at_psat(p_sat=p_ref) # [kg/m^3] Saturation density
    # Prandtl number
    Pr = fp.get_Prandtl(T=T, p=p_ref) # [-] Prandtl number
    Pr[0] = fp.get_saturation_Prandtl_gas(p_sat=p_ref) # [-] Saturation Prandtl
    # Thermal conductivity 
    kappa = fp.get_thermal_conductivity(T=T, p=p_ref) # [W/(m*K)] Conductivity
    kappa[0] = fp.get_gas_saturation_conductivity(p_sat=p_ref) # [W/(m*K)] Saturation conductivity
    # Viscosity
    mu = fp.get_viscosity(T=T, p=p_ref) # [Pa*s] Viscosity
    mu[0] = fp.get_gas_saturation_viscosity(p_sat=p_ref) # [Pa*s] Saturation viscosity
    return {\
        "T":T, # [K]
        "dT": dT, # [K]
        "rho": rho, # [kg/m^3]
        "h": h, # [J/kg]
        "Q_dot": Q_dot, # [W]
        "Pr": Pr, # [-]
        "kappa": kappa, # [W/(m*K)]
        "mu": mu, # [Pa*s]
        }
    
def full_homogenous_preparation(T_inlet, T_outlet, m_dot, p_ref, steps_l, steps_tp, steps_g, fp:FluidProperties):

    p_l = prepare_single_phase_liquid(
        T_inlet=T_inlet,
        steps=steps_l,
        p_ref=p_ref,
        m_dot=m_dot,
        fp=fp
    )

    p_tp = prepare_homogenous_transition(
        p=p_ref,
        m_dot=m_dot,
        steps=steps_tp,
        fp = fp
    )

    p_g = prepare_single_phase_gas(
        T_outlet=T_outlet,
        steps=steps_g,
        p_ref=p_ref,
        m_dot=m_dot,
        fp=fp
    )  
    ## Dictionary of prepared values
    return {
        'p_l': p_l,
        'p_tp': p_tp,
        'p_g': p_g
        }

def full_homogenous_calculation(prepared_values, Nusselt_relations, A_channel, wetted_perimeter, D_hydraulic, m_dot, T_wall, p_ref, fp: FluidProperties, pressure_drop_relations=None, area_ratio_contraction=None):
    # Unpack prepared values
    p_l = prepared_values['p_l']
    p_tp = prepared_values['p_tp']
    p_g = prepared_values['p_g']

    res_l = calc_channel_single_phase(
        T = p_l['T'],
        Q_dot= p_l['Q_dot'],
        rho = p_l['rho'],
        Pr = p_l['Pr'],
        kappa = p_l['kappa'],
        mu = p_l['mu'],
        p_ref=p_ref,
        m_dot=m_dot,
        T_wall=T_wall,
        D_hydr=D_hydraulic,
        wetted_perimeter=wetted_perimeter,
        A_channel=A_channel,
        Nu_func=Nusselt_relations['Nu_func_liquid'],
        fp=fp
    )

    res_tp = calc_homogenous_transition(
        p_sat=p_ref,
        x=p_tp['x'],
        alpha=p_tp['alpha'],
        T_sat=p_tp['T_sat'],
        rho_l=p_tp['rho_l'],
        rho_g=p_tp['rho_g'],
        rho=p_tp['rho'],
        m_dot=m_dot,
        mu_l=p_tp['mu_l'],
        mu=p_tp['mu'],
        Pr_l=p_tp['Pr_l'],
        Pr=p_tp['Pr'],
        kappa_l=p_tp['kappa_l'],
        kappa=p_tp['kappa'],
        Q_dot=p_tp['Q_dot'],
        T_wall=T_wall,
        D_hydr=D_hydraulic,
        wetted_perimeter=wetted_perimeter,
        A_channel=A_channel,
        Nu_func_tp=Nusselt_relations['Nu_func_two_phase'],
        Nu_func_le=Nusselt_relations['Nu_func_le'],
        Nu_func_dryout=Nusselt_relations['Nu_func_dryout'],
        fp=fp
    )

    res_g = calc_channel_single_phase(
        T=p_g['T'],
        Q_dot=p_g['Q_dot'],
        rho=p_g['rho'],
        Pr=p_g['Pr'],
        kappa=p_g['kappa'],
        mu=p_g['mu'],
        p_ref=p_ref,
        m_dot=m_dot,
        T_wall=T_wall,
        D_hydr=D_hydraulic,
        wetted_perimeter=wetted_perimeter,
        A_channel=A_channel,
        Nu_func=Nusselt_relations['Nu_func_gas'],
        fp=fp
    )

    # Calculate the total length of all sections, which is the sum of the last elements in each length array
    L_total = res_l['L'][-1] + res_tp['L'][-1] + res_g['L'][-1] # [m]
    print("\n -- Lengths -- :")
    print("Liquid: {:4.3f} mm".format(res_l['L'][-1]*1e3))
    print("Two-phase: {:4.3f} mm".format(res_tp['L'][-1]*1e3))
    print("Gas: {:4.3f} mm".format(res_g['L'][-1]*1e3))
    # If necessary, calculate the pressure drop. If just ONE pressure drop relation is given, it stops being None.
    dP_total = None # [Pa] The pressure drop is returned as None if none is caluclated
    p_chamber = None # [Pa] The chamber pressure is returned as None if no pressure drop is calculated. It is purposefully not set to p_ref or p_inlet, as the user must be aware of which value he uses
    
    # NOTE: Contraction losses are calculated before the channels are determined, as it has quite an impact on heat transfer/channel length
    dP_contraction = None # [Pa] Pressure drop due to sudden contraction at start of channel
    dP_l = None # [Pa] Pressure drop in liquid section of channel
    dP_tp = None # [Pa] Pressure drop in two-phase section of channel
    dP_g = None # [Pa] Pressure drop in gas section of channel
    if (pressure_drop_relations is not None):
        # Check if contraction relation is given
        if 'contraction' in pressure_drop_relations:
            total_dynamic_pressure = 0.5*p_l['rho'][0]*res_l['u'][0]**2 # [Pa] Total dynamic pressure
            args_contraction= {
                'area_ratio_contraction': area_ratio_contraction,
                'Re_Dh_downstream': res_l['Re'][0], # [-] The downstream Reynold's number is simply the first element in the liquid phase Reynold's array
                'total_dynamic_pressure': total_dynamic_pressure , 
            }
            dP_contraction = pressure_drop_relations['contraction'](args=args_contraction) # [Pa] Pressure drop due to sudden contraction
        else:
            dP_contraction = 0 # [Pa] Pressure drop due to sudden contraction
        
        # Check if liquid relation is given
        if 'l' in pressure_drop_relations: # Calculate pressure drop if given
            args_l = {
                'Re' : res_l['Re'], # [-] Reynolds number for x=0
                'rho': p_l['rho'], # [kg/m^3] Density
                'u': res_l['u'], # [m/s] Velocity
                'D_hydraulic': D_hydraulic, # [m] Hydraulic diameter
                'delta_L': res_l['delta_L'], # Size of meshes
            }
            dP_l = pressure_drop_relations['l'](args=args_l) # [Pa] Pressure drop in liquid section of channel
        else: # No relation given, means it is assumed to be neglected
            dP_l = 0 # [Pa]

        # Check if two-phase relation is given
        if 'tp' in pressure_drop_relations:
        ## Two-phase-arguments
            args_tp = {
                'x': p_tp['x'], # [-] Array of vapour quality
                'Re_tp_l' : res_tp['Re'][0], # [-] Reynolds number for x=0
                'Re_tp_g' : res_tp['Re'][-1],   # [-] Similary, x=1
                'rho_tp_l': p_tp['rho'][0], # [kg/m^3] Liquid density when x=0
                'rho_tp_g': p_tp['rho'][-1], # [kg/m^3] Gas density when x=1
                'u_tp_l': res_tp['u'][0], # [m/s] Liquid Velocity for x=0
                'u_tp_g': res_tp['u'][-1], # [m/s] Gas Velocity for x=1
                'D_hydraulic': D_hydraulic, # [m] Hydraulic diameter
                'delta_L': res_tp['delta_L'], # Size of meshes
            }           
            dP_tp = pressure_drop_relations['tp'](args=args_tp) # [Pa] Pressure drop in two-phase seciton of channel
        else: # No relation given, means it is assumed to be neglected
            dP_tp = 0 # [Pa] Pressure drop in two-phase seciton of channel

        # Check if gas-relation is given

        if 'g' in pressure_drop_relations:
            # Gas arguments
            args_g = {
                'Re' : res_g['Re'], # [-] Reynolds number for x=0
                'rho': p_g['rho'], # [kg/m^3] Density
                'u': res_g['u'], # [m/s] Velocity
                'D_hydraulic': D_hydraulic, # [m] Hydraulic diameter
                'delta_L': res_g['delta_L'], # Size of meshes
            }
            # Calculate pressure drop in gas section 
            dP_g = pressure_drop_relations['g'](args=args_g)  # [Pa] Pressure drop in gas section of channel
        else:
            dP_g = 0 # [Pa] Pressure drop in gas section of channel


        # Summing up the results of calculated pressure drops in three sections of channel
        # IMPORTANT: Contraction loss is NOT added included in this total. It is already subtracted from p_chamber
        dP_total = dP_l + dP_tp + dP_g + dP_contraction # [Pa] Total pressure drop of whatever source
        p_chamber = p_ref - dP_total # [Pa] P_ref is assumed to be inlet pressure

        # Some values that may be of use, and are best calculated once correctly
        T_chamber = p_g['T'][-1] # Last temperature in prepared gas phase values is chamber temperature (according to IRT)


    return {
        'p_g': p_g,
        'p_tp': p_tp,
        'p_l': p_l,
        'res_l': res_l,
        'res_tp': res_tp,
        'res_g': res_g,
        'L_total': L_total,
        'dP_l': dP_l,
        'dP_tp': dP_tp,
        'dP_g': dP_g,
        'dP_contraction': dP_contraction,
        'dP_total': dP_total,
        'p_chamber': p_chamber,
        'T_chamber': T_chamber,
    }

def rectangular_multi_channel_homogenous_calculation(channel_amount, prepared_values, Nusselt_relations, pressure_drop_relations, w_channel, h_channel, m_dot, T_wall, p_inlet, fp: FluidProperties, area_ratio_contraction=None):
    """Same as full_homogenous_calculation, but with channels combined

    Args:
        channel_amount (-): Number of channels
        prepared_values (dict): Dictionary with prepared thermodynamic variables
        Nusselt_relations (dict): Nusselt relations for different phases
        pressure_drop_relations (dict): Pressure drop relations for different phases and effects
        A_channel (m^2): Area of a single channel
        wetted_perimeter (m): Circumference of a single channel (fot heat flux purposes)
        D_hydraulic (m): Reference for Reynolds numbers etc...
        m_dot (kg/s): mass flow
        T_wall (K): Wall temperature
        p_ref (Pa): Pressure through channel
        fp (FluidProperties): Object to access propellant properties with

    Returns:
        Results {}: Return results which are relevant for optimization
    """
    m_dot_channel = m_dot / channel_amount # [kg/s] Mass flow for a single channel
    A_channel = w_channel * h_channel # [m^2] Area of a single channel
    D_hydraulic = hydraulic_diameter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] 
    wetted_perimeter = wetted_perimeter_rectangular(w_channel=w_channel, h_channel=h_channel) # [m] 

    # dP_contraction = None # [Pa] Contraction pressure drop is set to None by default, so it raises obvious errors if one attemps to use it without calculating it before (0 would hide errors)
    # # Contraction losses must be calculated here, as this pressure drop does significantly affect heat transfer
    # if pressure_drop_relations is not None:
    #     if 'contraction'  in pressure_drop_relations:
    #         area_ratio_contraction = w_channel*channel_amount / w_inlet_manifold # [-]
    #         dP_contraction = calc_contraction_loss(dP_func=pressure_drop_relations['contraction'])
            

    

    # Simulation of the flow resulting in desired channel length
    res = full_homogenous_calculation(
            prepared_values=prepared_values,
            Nusselt_relations=Nusselt_relations,
            pressure_drop_relations=pressure_drop_relations,
            A_channel=A_channel,
            wetted_perimeter=wetted_perimeter,
            D_hydraulic=D_hydraulic,
            m_dot=m_dot_channel,
            T_wall=T_wall,
            p_ref=p_inlet,
            fp=fp,
            area_ratio_contraction=area_ratio_contraction
            )

    return res


def calc_two_phase_frictional_pressure_drop_low_Reynolds(args):
    
    # Müller-Steinhagen and Heck approach
    # Pressure drop per meter if fully liquid flow (x=0)
    A = 64/args['Re_tp_l'] * 1/(2*args['D_hydraulic']) * args['rho_tp_l'] * args['u_tp_l']**2 # [Pa/m] Pressure drop per unit length if fully liquid
    # Pressure drop per meter if fully gaseous
    B = 64/args['Re_tp_g'] * 1/(2*args['D_hydraulic']) * args['rho_tp_g'] * args['u_tp_g']**2 # [Pa/m] pressure drop per unit length if fully gaseous

    dp_dl = ( A + 2 * (B-A) * args['x'] )*(1-args['x'])**(1/3) + B*args['x']**3 # [Pa/m] Pressure drop for two-phase flow
    # Each section has varying lengths, so remember dL is probably not constant
    return np.sum(dp_dl * args['delta_L']) # [Pa] Two-phase frictional pressure drop

def calc_single_phase_frictional_pressure_drop_low_Reynolds(args):
    # Müller-Steinhagen and Heck approach for Reynolds below Re << 1187
    dp_dl= 64/args['Re'] * 1/(2*args['D_hydraulic']) * args['rho'] * args['u']**2 # [Pa/m] Pressure drop per unit length if fully liquid
    return np.sum(dp_dl * args['delta_L']) # [Pa] Two-phase frictional pressure drop

def calc_single_phase_frictional_pressure_drop_high_Reynolds(args):
    # Müller-Steinhagen and Heck approach for Reynolds above Re > 1187
    dp_dl= 0.3164/args['Re']**(1/4) * 1/(2*args['D_hydraulic']) * args['rho'] * args['u']**2 # [Pa/m] Pressure drop per unit length if fully liquid
    return np.sum(dp_dl * args['delta_L']) # [Pa] Two-phase frictional pressure drop

def calc_single_phase_contraction_pressure_drop_Kawahara2015(args):
    # Kawahara2015 contraction pressure drop for Reynold be Re < 2200

    B = 0.0645 * np.log(args['Re_Dh_downstream']) - 0.00792
    C_C = B + ( 1 - B)*args['area_ratio_contraction']**(4.519)
    zeta_C = (1-1/C_C)**2
    
    return zeta_C * args['total_dynamic_pressure']