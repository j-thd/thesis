""" In this file the main optimization algorithm is implemented, and it takes the following steps:

    * Calculate chamber parameters from the defined top-level-requirements and fixed parameters of the design space


"""

import numpy as np
from scipy.optimize import root_scalar

import basic.chamber as chamber
from thermo.prop import FluidProperties
from basic.IRT import engine_performance_from_F_and_T
import models.one_D as oneD
from models.one_D import BottomWallTemperatureTooLowError 

def optim_P_total(channel_amount, w_channel, h_channel, inlet_manifold_length_factor, inlet_manifold_width_factor, l_exit_manifold, w_channel_spacing, w_outer_margin, T_wall, p_ref, m_dot, \
    prepared_values, Nusselt_relations, pressure_drop_relations, convergent_half_angle, divergent_half_angle, F_desired, p_back, AR_exit, emissivity_top, fp: FluidProperties, wall_args=None):
    # First calculate area ratio at entrance so it can be provided to contraction pressure drop function
    w_inlet_manifold = chamber.inlet_manifold_width(\
            channel_amount=channel_amount,
            w_channel=w_channel,
            w_channel_spacing=w_channel_spacing,
            inlet_manifold_width_factor=inlet_manifold_width_factor) # [m] Total width of inlet manifold
    area_ratio_contraction = chamber.area_ratio_contraction(\
        w_inlet_manifold=w_inlet_manifold,
        w_channel=w_channel,
        channel_amount=channel_amount) # [-] Area ratio of sudden contraction
    


    def x(T):
        wall_args['T_wall_bottom'] = T # [K] Add the bottom wall temperature in the wall argument dictionary
        # Calculate heat transfer to determine channel length
        res = oneD.rectangular_multi_channel_homogenous_calculation(\
            channel_amount=channel_amount,
            prepared_values=prepared_values,
            Nusselt_relations=Nusselt_relations,
            pressure_drop_relations=pressure_drop_relations,
            w_channel=w_channel,
            h_channel=h_channel,
            m_dot=m_dot,
            T_wall=T_wall,
            p_inlet=p_ref,
            fp=fp,
            area_ratio_contraction=area_ratio_contraction,
            wall_args=wall_args)
        
        #print(res['Q_total_bottom_plane_heat_balance'])
        return res['Q_total_bottom_plane_heat_balance']

    ## ROOT FINDING PROBLEM FOR T_WALL_BOTTOM
    # The proper bounds for the optimization are very sensitive to wall temperature and chamber temperature
    # It is hard to pick bounds that are valid for a wide range of those parameters
    # The problem is that that a too low lower bound results in model results that are wholly invalid and break the optimization
    # Too high a lower bound results in the solution not being within the bounds at all.
    # The quick and dirty solution here is to stepwise move the bounds down until the valid solution is between it.
    # The stepsize must be small enough to reliably avoid the scenario where the calculations are completely incorrect.
    delta_T_bounds = (T_wall+prepared_values['p_l']['T'][0])/2 # [K] The step size is half the top wall temperature and the inlet temperature initially
   
    ## The resultant channel length, pressure drops, etc, depends on bottom wall temperature. It must be iterated towards here.
    T_wall_bottom_min = T_wall-delta_T_bounds# [K] Minimum temperature is in practice probably not lower than chamber temperature (not impossible, though)
    T_wall_bottom_max = T_wall
    bounds_are_correct = False # [bool] Is False until a solution has been found
    # Setting the proper bounds is already hard, and here, bifurcation will be used.
    # The bounds are correct once they get the opposite sign
    ## NOTE: this while statement is a bit too complex and improper, but was a late and necessary fix to the code
    while not bounds_are_correct:
        # If one of the temperatures is too low, then the bounds must be raised.
        try:
            a = x(T_wall_bottom_max) 
            b = x(T_wall_bottom_min)
        except BottomWallTemperatureTooLowError:
            # This means the lower bound b is too low. It is raised by half the previous step size
            delta_T_bounds = delta_T_bounds / 2
            T_wall_bottom_min += delta_T_bounds
            continue # The loop must be aborted if this error was raised.
        
        # If the solutions have the same sign:
        if np.sign(a) == np.sign(b):
            T_wall_bottom_max = T_wall_bottom_min # Maker the upper bound the previous lower bound.
            # Half the interval again by halving the delta T
            delta_T_bounds = delta_T_bounds / 2
            T_wall_bottom_min -= delta_T_bounds 
            
        else:
            bounds_are_correct = True

    root_result = root_scalar(x, bracket=([T_wall_bottom_min,T_wall_bottom_max]))
            
        

    if root_result.converged:
        T_wall_bottom = root_result.root # [m^2]
        wall_args['T_wall_bottom'] = T_wall_bottom
        #print(" T_wall_bottom: {:3.2f} K".format(T_wall_bottom))
    else:
        raise Exception("No solution found for given temperature")
    
    res = oneD.rectangular_multi_channel_homogenous_calculation(\
            channel_amount=channel_amount,
            prepared_values=prepared_values,
            Nusselt_relations=Nusselt_relations,
            pressure_drop_relations=pressure_drop_relations,
            w_channel=w_channel,
            h_channel=h_channel,
            m_dot=m_dot,
            T_wall=T_wall,
            p_inlet=p_ref,
            fp=fp,
            area_ratio_contraction=area_ratio_contraction,
            wall_args=wall_args)

    # Check if these solutions are eventually valid
    assert(np.all(res['res_l']['delta_L']>=0))
    assert(np.all(res['res_tp']['delta_L']>=0))
    assert(np.all(res['res_g']['delta_L']>=0))
        
    l_channel = res['L_total'] # [m] Channel length
    p_chamber = res['p_chamber'] # [Pa] Pressure out channel outlet/inlet of nozzle according to IRT
    T_chamber = res['T_chamber']

    # Nozzle performance must be recalculated, if pressure drop occured
    assert(pressure_drop_relations is not None) # Just check is one was calculated

    if p_chamber <= 0: # If pressure drop was so large p_chamber was negative, the solution is invalid, and NaNs must be returned
        # Return the same dictionary, but all NaN.
        return {
            'P_loss': np.nan, # [W] Current approximation of heat loss
            'l_channel': np.nan, # [m] Channel length
            'pressure_drop': np.nan, # [Pa] Total pressure drop
            'dP_contraction': np.nan, # [Pa] Pressure drop de to contraction
            'w_total': np.nan, # [m] Total chip width
            'l_total': np.nan, # [m] Total chip length
            'is_channel_choked': np.nan,
            'Re_channel_exit': np.nan, # [-] Reynolds number at channel exit
            'pressure_drop_punishment': 0-p_chamber# [-] This factor is used to help avoid the optimization converging at invalid solutions with a constant punishment factor
    }
    # This only runs if pressure drop was not too large.
    # In that case the throat area must be recalculated still
    else:
            

        ep = engine_performance_from_F_and_T(
                F_desired=F_desired,
                p_chamber=p_chamber,
                T_chamber=T_chamber,
                AR_exit=AR_exit,
                p_back=p_back,
                fp=fp
            )

        # New throat widtl_inlet_manifoldmined from desired performance parameter
        w_throat_new = ep['A_throat']/h_channel # [m]
        # Check if new pressure, does not cause choked heating channels
        is_channel_choked = (w_channel*h_channel*channel_amount < ep['A_throat']) # [bool] If combined channel area is smaller than throat, they are choked.
        # Check new Reynolds number at channel exit for laminar/turbulent flow assumptions
        hydraulic_diameter = chamber.hydraulic_diameter_rectangular(w_channel=w_channel,h_channel=h_channel)
        A_channel = h_channel*w_channel # [m^2] Cross-sectional area of channel
        Re_channel_l = res['res_l']['Re'][0] 
        Re_channel_tp = res['res_tp']['Re'][0]
        Re_channel_g = res['res_g']['Re'][0]
        M_channel_exit_after_dP = fp.get_Mach_number(T=T_chamber, p=p_chamber, u=res['res_g']['u'][-1])

        Re_channel_exit_after_dP = fp.get_Reynolds_from_mass_flow(T=T_chamber, p=p_chamber,L_ref=hydraulic_diameter, m_dot=m_dot, A=A_channel)

        l_outlet = chamber.outlet_length(\
            w_channel=w_channel,
            w_channel_spacing=w_channel_spacing,
            channel_amount=channel_amount,
            convergent_half_angle=convergent_half_angle,
            w_throat=w_throat_new,
            divergent_half_angle=divergent_half_angle,
            AR_exit=AR_exit,
            l_exit_manifold=l_exit_manifold
            )
        
        
        l_inlet_manifold = chamber.inlet_manifold_length(w_inlet_manifold=w_inlet_manifold, inlet_manifold_length_factor=inlet_manifold_length_factor)
        l_total = chamber.total_chip_length(l_inlet_manifold=l_inlet_manifold, l_channel=l_channel, l_outlet=l_outlet)
        w_total = chamber.total_chip_width(\
            w_inlet_manifold=w_inlet_manifold,
            w_outer_margin=w_outer_margin,
            w_throat=w_throat_new,
            AR_exit=AR_exit)

        A_chip = l_total*w_total # [m^2]

        P_rad_loss_top = chamber.radiation_loss(T=T_wall, A=A_chip, emmisivity=emissivity_top) # [W]  Radiation loss through top of chip
        P_rad_loss_bottom = chamber.radiation_loss(T=wall_args['T_wall_bottom'],A=A_chip, emmisivity=wall_args['emissivity_chip_bottom']) # Radiation loss through bottom of chip
        P_loss = P_rad_loss_top + P_rad_loss_bottom # [W] Total heat loss
        return {'P_loss': P_loss, # [W] Current approximation of heat loss
                'l_channel': l_channel, # [m] Channel length
                'pressure_drop': res['dP_total'], # [Pa] Total pressure drop
                'dP_contraction': res['dP_contraction'], # [Pa] Pressure drop due to contraction
                'w_total': w_total, # [m] Total chip width
                'l_total': l_total, # [m] Total chip length
                'is_channel_choked': is_channel_choked, # [bool] 0 and 1 in this case, to check if channels are not too small
                'Re_channel_exit_after_dP': Re_channel_exit_after_dP, # [-] Check Reynolds number to see if laminar flow assumptions still hold
                'M_channel_exit_after_dP': M_channel_exit_after_dP, # [-] Check Mach number to see if it is a problem
                'Re_channel_l': Re_channel_l,
                'Re_channel_tp': Re_channel_tp,
                'Re_channel_g': Re_channel_g,
                'T_wall_bottom': T_wall_bottom, # [K] Bottom wall temperature 
        }